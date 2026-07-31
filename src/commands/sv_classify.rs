use crate::alignment_record::Strand;
use crate::impg::CigarOp;
use crate::impg_index::ImpgIndex;
use log::warn;
use rustc_hash::FxHashSet;
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
    pub tandem_cv_threshold: f32,
    pub merge_gap: u32,
    pub min_support: u32,
}

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

            // CIGAR-based gap extraction (DEL/INS/TDUP/TCON)
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

        // D events → DEL (consistent gap) or TCON (variable gap = tandem contraction)
        if !del_sizes.is_empty() {
            let support = del_sizes.len() as u32;
            let size = median(&del_sizes);
            if coefficient_of_variation(&del_sizes) > filters.tandem_cv_threshold {
                if size >= filters.tcon_min && size <= filters.tcon_max {
                    calls.push(SvCall {
                        target_name: target_name.to_string(),
                        target_start: locus_start,
                        target_end: locus_end_pos,
                        sv_type: SvType::Tcon,
                        size,
                        support,
                        query_regions: del_query_regions,
                    });
                }
            } else if size >= filters.del_min && size <= filters.del_max {
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

        // I events → INS (consistent) or TDUP (variable = tandem duplication)
        if !ins_sizes.is_empty() {
            let support = ins_sizes.len() as u32;
            let size = median(&ins_sizes);
            if coefficient_of_variation(&ins_sizes) > filters.tandem_cv_threshold {
                if size >= filters.tdup_min && size <= filters.tdup_max {
                    calls.push(SvCall {
                        target_name: target_name.to_string(),
                        target_start: locus_start,
                        target_end: locus_end_pos,
                        sv_type: SvType::Tdup,
                        size,
                        support,
                        query_regions: ins_query_regions,
                    });
                }
            } else if size >= filters.ins_min && size <= filters.ins_max {
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

fn coefficient_of_variation(values: &[u32]) -> f32 {
    if values.len() < 2 {
        return 0.0;
    }
    let n = values.len() as f32;
    let mean = values.iter().sum::<u32>() as f32 / n;
    if mean == 0.0 {
        return 0.0;
    }
    let var = values
        .iter()
        .map(|&v| {
            let d = v as f32 - mean;
            d * d
        })
        .sum::<f32>()
        / n;
    var.sqrt() / mean
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
