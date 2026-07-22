use crate::impg::CigarOp;
use crate::impg_index::ImpgIndex;
use std::io::{self, Write};

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
    pub vcf: bool,
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
}

struct GapEvent {
    target_start: i32,
    target_end: i32,
    size: u32,
    is_del: bool,
}

pub fn run(
    impg: &impl ImpgIndex,
    scan_regions: Vec<(u32, i32, i32)>,
    filters: &SvFilters,
) -> io::Result<()> {
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    if filters.vcf {
        write_vcf_header(&mut out)?;
    } else {
        writeln!(out, "#chrom\tstart\tend\tsv_type\tsize\tsupport")?;
    }

    for (target_id, t_start, t_end) in scan_regions {
        let target_name = impg.seq_index().get_name(target_id).unwrap().to_string();
        let results = impg.query(target_id, t_start, t_end, true, None, None, false)?;

        let mut gap_events: Vec<GapEvent> = Vec::new();
        let mut direct_calls: Vec<SvCall> = Vec::new();

        // result[0] is always the identity mapping of the target region itself
        for (query_iv, cigar_ops, target_iv) in results.iter().skip(1) {
            let query_id = query_iv.metadata;
            let target_span = (target_iv.last - target_iv.first) as u32;

            // INV: reverse-strand alignment encodes as query_iv.first > query_iv.last
            if query_iv.first > query_iv.last
                && target_span >= filters.inv_min
                && target_span <= filters.inv_max
            {
                direct_calls.push(SvCall {
                    target_name: target_name.clone(),
                    target_start: target_iv.first,
                    target_end: target_iv.last,
                    sv_type: SvType::Inv,
                    size: target_span,
                    support: 1,
                });
            }

            // TRA: alignment from a different chromosome (PanSN-aware)
            let query_name = impg.seq_index().get_name(query_id).unwrap();
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
                });
            }

            // CIGAR-based gap extraction (DEL/INS/TDUP/TCON)
            if !cigar_ops.is_empty() {
                gap_events.extend(extract_gap_events(cigar_ops, target_iv.first));
            }
        }

        for call in classify_gap_loci(&target_name, gap_events, filters) {
            if call.support >= filters.min_support {
                emit_call(&mut out, &call, filters.vcf)?;
            }
        }
        for call in &direct_calls {
            if call.support >= filters.min_support {
                emit_call(&mut out, call, filters.vcf)?;
            }
        }
    }

    Ok(())
}

fn extract_gap_events(cigar_ops: &[CigarOp], target_start: i32) -> Vec<GapEvent> {
    let mut events = Vec::new();
    let mut target_pos = target_start;

    for op in cigar_ops {
        match op.op() {
            'D' => {
                events.push(GapEvent {
                    target_start: target_pos,
                    target_end: target_pos + op.len(),
                    size: op.len() as u32,
                    is_del: true,
                });
                target_pos += op.len();
            }
            'I' => {
                events.push(GapEvent {
                    target_start: target_pos,
                    target_end: target_pos,
                    size: op.len() as u32,
                    is_del: false,
                });
                // I ops consume no target bases
            }
            _ => {
                target_pos += op.target_delta();
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

    for event in events {
        if !current.is_empty() && event.target_start > locus_end + filters.merge_gap as i32 {
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

fn write_vcf_header(out: &mut impl Write) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(
        out,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length in bp\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of supporting alignments\">"
    )?;
    writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
}

fn emit_call(out: &mut impl Write, call: &SvCall, vcf: bool) -> io::Result<()> {
    if vcf {
        writeln!(
            out,
            "{}\t{}\t.\tN\t<{}>\t.\tPASS\tSVTYPE={};SVLEN={};SUPPORT={}",
            call.target_name,
            call.target_start + 1, // VCF is 1-based
            call.sv_type.label(),
            call.sv_type.label(),
            call.size,
            call.support,
        )
    } else {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            call.target_name,
            call.target_start,
            call.target_end,
            call.sv_type.label(),
            call.size,
            call.support,
        )
    }
}
