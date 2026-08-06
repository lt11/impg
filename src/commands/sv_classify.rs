use crate::alignment_record::Strand;
use crate::impg::CigarOp;
use crate::impg_index::ImpgIndex;
use log::warn;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{self, Write};

/// Controls which query sequences are considered when calling variants
/// against a target locus.
pub enum QueryFilter {
    /// Default: use every query except ones belonging to the target
    ExcludeSameSample,
    /// Explicit set of allowed query sequence ids (from `--query-name`)
    Explicit(FxHashSet<u32>),
}

impl QueryFilter {
    fn allows(&self, query_id: u32, query_name: &str, target_name: &str) -> bool {
        match self {
            QueryFilter::ExcludeSameSample => sample_part(query_name) != sample_part(target_name),
            QueryFilter::Explicit(ids) => ids.contains(&query_id),
        }
    }
}

pub struct SvFilters {
    pub del_min: u32,
    pub del_max: u32,
    pub ins_min: u32,
    pub ins_max: u32,
    pub inv_min: u32,
    pub inv_max: u32,
    pub tra_min: u32,
    pub tra_max: u32,
    pub tdup_min: u32,
    pub tdup_max: u32,
    pub tcon_min: u32,
    pub tcon_max: u32,
    pub merge_gap: u32,
    pub min_support: u32,
    pub inv_proxy_tolerance_pct: u32,
    pub inv_exclusion_buffer: u32,
}

/// MUM&Co uses fixed 50 bp thresholds throughout its overlap-based
/// duplication/contraction detection: the minimum block overlap that counts
/// as evidence, and the maximum gap on the other axis for the pairing to be
/// considered "clean" (i.e. truly tandem rather than coincidental).
const BLOCK_MIN_OVERLAP: i32 = 50;
const BLOCKS_MAX_CLEAN_GAP: i32 = 50;

enum SvType {
    Del,
    Ins,
    Inv,
    Tra,
    Tdup,
    Tcon,
}

impl SvType {
    fn label(&self) -> &'static str {
        match self {
            SvType::Del => "DEL",
            SvType::Ins => "INS",
            SvType::Inv => "INV",
            SvType::Tra => "TRA",
            SvType::Tdup => "TDUP",
            SvType::Tcon => "TCON",
        }
    }
}

struct SvCall {
    target_name: String,
    target_start: i32,
    target_end: i32,
    sv_type: SvType,
    size: u32,
    support: u32,
    /// (query_name, query_start, query_end) for each alignment supporting this call
    query_regions: Vec<(String, i32, i32)>,
}

struct GapEvent {
    target_start: i32,
    target_end: i32,
    size: u32,
    is_del: bool,
    query_name: String,
    query_start: i32,
    query_end: i32,
}

/// One alignment chain (AdjustedInterval) between a query and the target
/// locus, used for MUM&Co-style overlap detection (TDUP/TCON). Unlike
/// `GapEvent`, this captures the whole block, not individual CIGAR ops.
struct AlignmentBlock {
    query_name: String,
    target_start: i32,
    target_end: i32,
    query_start: i32,
    query_end: i32,
    reverse: bool,
}

pub fn run(
    impg: &impl ImpgIndex,
    scan_regions: Vec<(u32, i32, i32)>,
    query_filter: &QueryFilter,
    filters: &SvFilters,
) -> io::Result<()> {
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    writeln!(
        out,
        "chrom\tstart\tend\tsv_type\tsize\tsupport\tquery_chrom\tquery_start\tquery_end"
    )?;

    for (target_id, t_start, t_end) in scan_regions {
        let target_name = impg
            .seq_index()
            .get_name(target_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("No sequence name found for target id {}", target_id),
                )
            })?
            .to_string();
        let results = impg.query(target_id, t_start, t_end, true, None, None, false)?;

        let mut gap_events: Vec<GapEvent> = Vec::new();
        let mut direct_calls: Vec<SvCall> = Vec::new();
        let mut blocks_by_query: FxHashMap<u32, Vec<AlignmentBlock>> = FxHashMap::default();

        // result[0] is always the identity mapping of the target region itself
        for (query_iv, cigar_ops, target_iv) in results.iter().skip(1) {
            let query_id = query_iv.metadata;
            let query_name = impg.seq_index().get_name(query_id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("No sequence name found for query id {}", query_id),
                )
            })?;

            if !query_filter.allows(query_id, query_name, &target_name) {
                continue;
            }

            let target_span = (target_iv.last - target_iv.first) as u32;
            let reverse = query_iv.first > query_iv.last;
            let (query_start, query_end) = if reverse {
                (query_iv.last, query_iv.first)
            } else {
                (query_iv.first, query_iv.last)
            };

            // INV: reverse-strand alignment encodes as query_iv.first > query_iv.last
            if reverse && target_span >= filters.inv_min && target_span <= filters.inv_max {
                direct_calls.push(SvCall {
                    target_name: target_name.clone(),
                    target_start: target_iv.first,
                    target_end: target_iv.last,
                    sv_type: SvType::Inv,
                    size: target_span,
                    support: 1,
                    query_regions: vec![(query_name.to_string(), query_start, query_end)],
                });
            }

            // TRA: alignment from a different chromosome (PanSN-aware)
            if chrom_part(&target_name) != chrom_part(query_name)
                && target_span >= filters.tra_min
                && target_span <= filters.tra_max
            {
                direct_calls.push(SvCall {
                    target_name: target_name.clone(),
                    target_start: target_iv.first,
                    target_end: target_iv.last,
                    sv_type: SvType::Tra,
                    size: target_span,
                    support: 1,
                    query_regions: vec![(query_name.to_string(), query_start, query_end)],
                });
            }

            // CIGAR-based gap extraction (DEL/INS)
            if !cigar_ops.is_empty() {
                let strand = if reverse {
                    Strand::Reverse
                } else {
                    Strand::Forward
                };
                gap_events.extend(extract_gap_events(
                    cigar_ops,
                    target_iv.first,
                    query_iv.first,
                    strand,
                    query_name,
                ));
            }

            // Whole-block record for MUM&Co-style TDUP/TCON overlap detection
            blocks_by_query
                .entry(query_id)
                .or_default()
                .push(AlignmentBlock {
                    query_name: query_name.to_string(),
                    target_start: target_iv.first,
                    target_end: target_iv.last,
                    query_start,
                    query_end,
                    reverse,
                });
        }

        // A query's own CIGAR frequently carries small compensating indels
        // right at (or inside) an inversion breakpoint it also produced —
        // alignment-seam noise from stitching the reverse-strand block back
        // onto its forward-oriented flank, not an independent SV. Drop that
        // query's DEL/INS gap events wherever they fall within its own
        // strand-based INV call's target span (plus a small buffer).
        let mut inv_spans_by_query: FxHashMap<&str, Vec<(i32, i32)>> = FxHashMap::default();
        for call in &direct_calls {
            if matches!(call.sv_type, SvType::Inv) {
                for (query_name, _, _) in &call.query_regions {
                    inv_spans_by_query
                        .entry(query_name.as_str())
                        .or_default()
                        .push((call.target_start, call.target_end));
                }
            }
        }
        if !inv_spans_by_query.is_empty() {
            let buffer = filters.inv_exclusion_buffer as i64;
            gap_events.retain(|e| {
                !inv_spans_by_query
                    .get(e.query_name.as_str())
                    .is_some_and(|spans| {
                        spans.iter().any(|&(s, t)| {
                            e.target_start as i64 >= s as i64 - buffer
                                && e.target_end as i64 <= t as i64 + buffer
                        })
                    })
            });
        }

        for call in classify_gap_loci(&target_name, gap_events, filters) {
            if call.support >= filters.min_support {
                emit_call(&mut out, &call)?;
            }
        }
        for call in &direct_calls {
            if call.support >= filters.min_support {
                emit_call(&mut out, call)?;
            }
        }
        // TDUP/TCON: for each query present at this locus, look for that
        // query's own consecutive alignment blocks overlapping each other
        // (MUM&Co's signal). Scoped per query to avoid flagging the routine
        // case of many unrelated queries all covering the same target span.
        for blocks in blocks_by_query.values() {
            for call in detect_tandem_dups(&target_name, blocks, filters) {
                if call.support >= filters.min_support {
                    emit_call(&mut out, &call)?;
                }
            }
            for call in detect_tandem_contractions(&target_name, blocks, filters) {
                if call.support >= filters.min_support {
                    emit_call(&mut out, &call)?;
                }
            }
        }
    }

    Ok(())
}

fn extract_gap_events(
    cigar_ops: &[CigarOp],
    target_start: i32,
    query_start: i32,
    strand: Strand,
    query_name: &str,
) -> Vec<GapEvent> {
    let mut events = Vec::new();
    let mut target_pos = target_start;
    let mut query_pos = query_start;

    for op in cigar_ops {
        match op.op() {
            'D' => {
                events.push(GapEvent {
                    target_start: target_pos,
                    target_end: target_pos + op.len(),
                    size: op.len() as u32,
                    is_del: true,
                    query_name: query_name.to_string(),
                    query_start: query_pos,
                    query_end: query_pos,
                });
                target_pos += op.len();
            }
            'I' => {
                // I ops consume no target bases, only query bases
                let before = query_pos;
                query_pos += op.query_delta(strand);
                let (q_start, q_end) = if before <= query_pos {
                    (before, query_pos)
                } else {
                    (query_pos, before)
                };
                events.push(GapEvent {
                    target_start: target_pos,
                    target_end: target_pos,
                    size: op.len() as u32,
                    is_del: false,
                    query_name: query_name.to_string(),
                    query_start: q_start,
                    query_end: q_end,
                });
            }
            _ => {
                target_pos += op.target_delta();
                query_pos += op.query_delta(strand);
            }
        }
    }

    events
}

fn classify_gap_loci(
    target_name: &str,
    mut events: Vec<GapEvent>,
    filters: &SvFilters,
) -> Vec<SvCall> {
    if events.is_empty() {
        return Vec::new();
    }

    events.sort_by_key(|e| e.target_start);

    // Sweep-merge events within merge_gap bp into per-locus groups
    let mut loci: Vec<Vec<GapEvent>> = Vec::new();
    let mut current: Vec<GapEvent> = Vec::new();
    let mut locus_end = i32::MIN;
    let mut warned_merge_overflow = false;

    for event in events {
        let merge_threshold = match locus_end.checked_add(filters.merge_gap as i32) {
            Some(threshold) => threshold,
            None => {
                if !warned_merge_overflow {
                    warn!(
                        "sv-classify: target '{target_name}' locus end \
                         ({locus_end}) + --merge-gap ({}) overflows i32; \
                         capping merge distance at i32::MAX for the rest \
                         of this target (results may over-merge)",
                        filters.merge_gap
                    );
                    warned_merge_overflow = true;
                }
                i32::MAX
            }
        };
        if !current.is_empty() && event.target_start > merge_threshold {
            loci.push(std::mem::take(&mut current));
            locus_end = i32::MIN;
        }
        locus_end = locus_end.max(event.target_end);
        current.push(event);
    }
    if !current.is_empty() {
        loci.push(current);
    }

    let mut calls = Vec::new();

    for locus in loci {
        let (inv_calls, locus) = extract_inv_proxies(target_name, locus, filters);
        calls.extend(inv_calls);

        if locus.is_empty() {
            continue;
        }

        let locus_start = locus.iter().map(|e| e.target_start).min().unwrap();
        let locus_end_pos = locus.iter().map(|e| e.target_end).max().unwrap();

        let del_sizes: Vec<u32> = locus.iter().filter(|e| e.is_del).map(|e| e.size).collect();
        let ins_sizes: Vec<u32> = locus.iter().filter(|e| !e.is_del).map(|e| e.size).collect();
        let del_query_regions: Vec<(String, i32, i32)> = locus
            .iter()
            .filter(|e| e.is_del)
            .map(|e| (e.query_name.clone(), e.query_start, e.query_end))
            .collect();
        let ins_query_regions: Vec<(String, i32, i32)> = locus
            .iter()
            .filter(|e| !e.is_del)
            .map(|e| (e.query_name.clone(), e.query_start, e.query_end))
            .collect();

        // D events → DEL
        if !del_sizes.is_empty() {
            let support = del_sizes.len() as u32;
            let size = median(&del_sizes);
            if size >= filters.del_min && size <= filters.del_max {
                calls.push(SvCall {
                    target_name: target_name.to_string(),
                    target_start: locus_start,
                    target_end: locus_end_pos,
                    sv_type: SvType::Del,
                    size,
                    support,
                    query_regions: del_query_regions,
                });
            }
        }

        // I events → INS
        if !ins_sizes.is_empty() {
            let support = ins_sizes.len() as u32;
            let size = median(&ins_sizes);
            if size >= filters.ins_min && size <= filters.ins_max {
                calls.push(SvCall {
                    target_name: target_name.to_string(),
                    target_start: locus_start,
                    target_end: locus_end_pos,
                    sv_type: SvType::Ins,
                    size,
                    support,
                    query_regions: ins_query_regions,
                });
            }
        }
    }

    calls
}

/// Some aligners represent a short inversion not as a separate
/// reverse-strand alignment block but as a paired indel within an otherwise
/// forward-strand chain: an I op immediately followed by a D op of roughly
/// the same size (the inverted segment gets "inserted" in the query and the
/// original orientation's copy gets "deleted" from the target). This pulls
/// that signature out of a locus's gap events and reports it as INV instead
/// of the misleading co-located DEL+INS pair.
///
/// Pairing is done per query: a query contributes an INV proxy only if it
/// has exactly one D and one I event at this locus with sizes within
/// `filters.inv_proxy_tolerance_pct` of each other. Anything that doesn't
/// match this exact pattern (no pairing, more than one D or I from the same
/// query, size mismatch beyond tolerance) is left untouched for the normal
/// DEL/INS classification that follows.
///
/// Returns the extracted INV calls plus the remaining, unpaired events.
fn extract_inv_proxies(
    target_name: &str,
    locus: Vec<GapEvent>,
    filters: &SvFilters,
) -> (Vec<SvCall>, Vec<GapEvent>) {
    let mut by_query: FxHashMap<String, (Vec<GapEvent>, Vec<GapEvent>)> = FxHashMap::default();
    for event in locus {
        let entry = by_query.entry(event.query_name.clone()).or_default();
        if event.is_del {
            entry.0.push(event);
        } else {
            entry.1.push(event);
        }
    }

    let mut inv_calls = Vec::new();
    let mut remaining = Vec::new();

    for (_, (mut dels, mut inss)) in by_query {
        if dels.len() == 1 && inss.len() == 1 {
            let del = &dels[0];
            let ins = &inss[0];
            let max_size = del.size.max(ins.size) as i64;
            let diff = (del.size as i64 - ins.size as i64).abs();
            let within_tolerance =
                diff * 100 <= filters.inv_proxy_tolerance_pct as i64 * max_size;

            if within_tolerance && del.size >= filters.inv_min && del.size <= filters.inv_max {
                inv_calls.push(SvCall {
                    target_name: target_name.to_string(),
                    target_start: del.target_start,
                    target_end: del.target_end,
                    sv_type: SvType::Inv,
                    size: del.size,
                    support: 1,
                    query_regions: vec![(ins.query_name.clone(), ins.query_start, ins.query_end)],
                });
                continue;
            }
        }
        remaining.append(&mut dels);
        remaining.append(&mut inss);
    }

    (inv_calls, remaining)
}

/// One axis-agnostic view of an `AlignmentBlock` for overlap detection: the
/// "primary" axis is the one blocks are sorted and checked for overlap on,
/// the "secondary" axis is checked for a small ("clean") gap between the
/// two overlapping blocks. TDUP sorts on target/checks query; TCON sorts on
/// query/checks target — the geometry (including strand handling) is
/// otherwise identical, so both share this one implementation.
struct OverlapBlock {
    primary_start: i32,
    primary_end: i32,
    secondary_start: i32,
    secondary_end: i32,
    reverse: bool,
}

struct TandemEvent {
    /// Overlap region on the primary (sort) axis.
    primary_span: (i32, i32),
    /// Junction region on the secondary axis.
    secondary_span: (i32, i32),
    /// Magnitude of the primary-axis overlap.
    size: u32,
}

/// MUM&Co's overlap-based tandem detection: sort blocks by the primary axis,
/// and for each same-strand consecutive pair whose primary-axis ranges
/// overlap by at least `min_overlap`, check that the two blocks are also
/// closely adjacent (within `max_gap`) on the secondary axis — i.e. this
/// isn't just two unrelated blocks that happen to overlap on one axis, but
/// a genuinely tandem arrangement of the same underlying repeat.
fn detect_tandem_overlaps(
    blocks: &[OverlapBlock],
    min_overlap: i32,
    max_gap: i32,
) -> Vec<TandemEvent> {
    let mut sorted: Vec<&OverlapBlock> = blocks.iter().collect();
    sorted.sort_by_key(|b| b.primary_start);

    let mut events = Vec::new();
    for pair in sorted.windows(2) {
        let prev = pair[0];
        let curr = pair[1];
        if prev.reverse != curr.reverse {
            continue;
        }

        let overlap = prev.primary_end - curr.primary_start;
        if overlap < min_overlap {
            continue;
        }

        // Strand determines which ends of the secondary axis should be
        // adjacent: forward blocks correlate positively across the two
        // axes, reverse blocks correlate negatively.
        let gap = if prev.reverse {
            prev.secondary_start - curr.secondary_end
        } else {
            curr.secondary_start - prev.secondary_end
        };
        if gap.abs() > max_gap {
            continue;
        }

        let (sec_start, sec_end) = if prev.reverse {
            (
                curr.secondary_end.min(prev.secondary_start),
                curr.secondary_end.max(prev.secondary_start),
            )
        } else {
            (
                prev.secondary_end.min(curr.secondary_start),
                prev.secondary_end.max(curr.secondary_start),
            )
        };

        events.push(TandemEvent {
            primary_span: (curr.primary_start, prev.primary_end),
            secondary_span: (sec_start, sec_end),
            size: overlap as u32,
        });
    }
    events
}

/// TDUP: the same query's alignment blocks overlap in target coordinates —
/// the query carries an extra tandem copy of that target segment.
fn detect_tandem_dups(
    target_name: &str,
    blocks: &[AlignmentBlock],
    filters: &SvFilters,
) -> Vec<SvCall> {
    if blocks.len() < 2 {
        return Vec::new();
    }
    let query_name = &blocks[0].query_name;
    let overlap_blocks: Vec<OverlapBlock> = blocks
        .iter()
        .map(|b| OverlapBlock {
            primary_start: b.target_start,
            primary_end: b.target_end,
            secondary_start: b.query_start,
            secondary_end: b.query_end,
            reverse: b.reverse,
        })
        .collect();

    detect_tandem_overlaps(&overlap_blocks, BLOCK_MIN_OVERLAP, BLOCKS_MAX_CLEAN_GAP)
        .into_iter()
        .filter(|e| e.size >= filters.tdup_min && e.size <= filters.tdup_max)
        .map(|e| SvCall {
            target_name: target_name.to_string(),
            target_start: e.primary_span.0,
            target_end: e.primary_span.1,
            sv_type: SvType::Tdup,
            size: e.size,
            support: 1,
            query_regions: vec![(query_name.clone(), e.secondary_span.0, e.secondary_span.1)],
        })
        .collect()
}

/// TCON: the same query's alignment blocks overlap in query coordinates —
/// the target carries an extra tandem copy that this query is missing.
fn detect_tandem_contractions(
    target_name: &str,
    blocks: &[AlignmentBlock],
    filters: &SvFilters,
) -> Vec<SvCall> {
    if blocks.len() < 2 {
        return Vec::new();
    }
    let query_name = &blocks[0].query_name;
    let overlap_blocks: Vec<OverlapBlock> = blocks
        .iter()
        .map(|b| OverlapBlock {
            primary_start: b.query_start,
            primary_end: b.query_end,
            secondary_start: b.target_start,
            secondary_end: b.target_end,
            reverse: b.reverse,
        })
        .collect();

    detect_tandem_overlaps(&overlap_blocks, BLOCK_MIN_OVERLAP, BLOCKS_MAX_CLEAN_GAP)
        .into_iter()
        .filter(|e| e.size >= filters.tcon_min && e.size <= filters.tcon_max)
        .map(|e| SvCall {
            target_name: target_name.to_string(),
            target_start: e.secondary_span.0,
            target_end: e.secondary_span.1,
            sv_type: SvType::Tcon,
            size: e.size,
            support: 1,
            query_regions: vec![(query_name.clone(), e.primary_span.0, e.primary_span.1)],
        })
        .collect()
}

/// Extract the chromosome identifier from a sequence name.
/// For PanSN format (sample#haplotype#chromosome), returns the chromosome part.
/// For plain names, returns the full name.
fn chrom_part(name: &str) -> &str {
    match name.rfind('#') {
        Some(pos) => &name[pos + 1..],
        None => name,
    }
}

/// Extract the sample identifier from a sequence name.
/// For PanSN format (sample#haplotype#chromosome), returns the sample part.
/// For plain names, returns the full name.
fn sample_part(name: &str) -> &str {
    match name.find('#') {
        Some(pos) => &name[..pos],
        None => name,
    }
}

fn median(values: &[u32]) -> u32 {
    let mut s = values.to_vec();
    s.sort_unstable();
    let mid = s.len() / 2;
    if s.len() % 2 == 0 {
        (s[mid - 1] + s[mid]) / 2
    } else {
        s[mid]
    }
}

/// Emits one row per supporting query alignment (one-row-per-event), all
/// sharing the same locus-level chrom/start/end/sv_type/size/support.
fn emit_call(out: &mut impl Write, call: &SvCall) -> io::Result<()> {
    for (query_name, query_start, query_end) in &call.query_regions {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            call.target_name,
            call.target_start,
            call.target_end,
            call.sv_type.label(),
            call.size,
            call.support,
            query_name,
            query_start,
            query_end,
        )?;
    }
    Ok(())
}
