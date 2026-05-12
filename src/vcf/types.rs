//! VCF data types
//! This module use [ Variant Call Format(VCF) Version 4.2 Specification ](https://samtools.github.io/hts-specs/VCFv4.2.pdf) as standard reference

#[derive(Debug, Default)]
pub struct Contig {
    pub id: String,
    pub length: Option<u64>,
    pub assembly: Option<String>,
}

#[derive(Debug, Clone, Default)]
pub struct FormatDef {
    pub id: String,
    pub number: String,
    pub type_: String,
    pub description: String,
}

#[derive(Debug, Default)]
pub struct InfoDef {
    pub id: String,
    pub number: String,
    pub type_: String,
    pub description: String,
}

#[derive(Debug, Default)]
pub struct VcfMeta {
    pub contigs: Vec<Contig>,
    pub formats: Vec<FormatDef>,
    pub infos: Vec<InfoDef>,
    pub samples: Vec<String>,
}

impl Contig {
    pub fn parse(line: &str) -> Option<Self> {
        let content = line
            .strip_prefix("##contig=")?
            .strip_prefix('<')?
            .strip_suffix('>')?;
        let mut c = Self::default();
        for part in content.split(',') {
            let mut it = part.split('=');
            let key = it.next()?.trim();
            let val = it.next()?.trim();
            match key {
                "ID" => c.id = val.to_string(),
                "length" => c.length = val.parse().ok(),
                "assembly" => c.assembly = Some(val.to_string()),
                _ => {}
            }
        }
        Some(c)
    }
}

impl FormatDef {
    pub fn parse(line: &str) -> Option<Self> {
        let content = line
            .strip_prefix("##FORMAT=")?
            .strip_prefix('<')?
            .strip_suffix('>')?;
        let mut f = Self::default();
        for part in content.split(',') {
            let mut it = part.split('=');
            let key = it.next()?.trim();
            let val = it.next()?.trim();
            match key {
                "ID" => f.id = val.to_string(),
                "Number" => f.number = val.to_string(),
                "Type" => f.type_ = val.to_string(),
                "Description" => f.description = val.trim_matches('"').to_string(),
                _ => {}
            }
        }
        Some(f)
    }
}

use std::collections::HashMap;

use arrow::array::ArrayRef;

#[derive(Debug)]
pub struct VcfSample {
    pub format_def: FormatDef,
    pub array: ArrayRef,
}

/// VCF-to-Arrow conversion result struct
#[derive(Debug)]
pub struct VcfParseResult {
    pub meta: VcfMeta,

    /// CHROM - chromosome: An identifier from the reference genome or an angle-bracketed ID String ("\<ID\>")
    /// pointing to a contig in the assembly file (cf. the ##assembly line in the header). All entries for a specific
    /// CHROM should form a contiguous block within the VCF file.
    ///
    /// (String, no whitespace permitted, Required).
    pub chrom: ArrayRef,

    /// POS - position: The reference position, with the 1st base having position 1. Positions are sorted numerically,
    /// in increasing order, within each reference sequence CHROM. It is permitted to have multiple records with the same POS.
    /// Telomeres are indicated by using positions 0 or N+1, where N is the length of the correspondin chromosome or contig.
    ///
    /// (Integer, Required)
    pub pos: ArrayRef,

    /// ID - identifier: Semicolon-separated list of unique identifiers where available. If this is a dbSNP variant it is
    /// encouraged to use the rs number(s). No identifier should be present in more than one data record. If there is
    /// no identifier available, then the missing value should be used. (String, no whitespace or semicolons permitted)
    pub id: ArrayRef,

    /// REF - reference base(s): Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted.
    /// The value in the POS field refers to the position of the first base in the String. For simple insertions and
    /// deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT
    /// Strings must include the base before the event (which must be reflected in the POS field), unless the event
    /// occurs at position 1 on the contig in which case it must include the base after the event; this padding base is
    /// not required (although it is permitted) for e.g. complex substitutions or other events where all alleles have at
    /// least one base represented in their Strings. If any of the ALT alleles is a symbolic allele (an angle-bracketed
    /// ID String ("\<ID\>") then the padding base is required and POS denotes the coordinate of the base preceding
    /// the polymorphism. Tools processing VCF files are not required to preserve case in the allele Strings. (String,
    /// Required).
    pub _ref: ArrayRef,

    /// ALT - alternate base(s): Comma separated list of alternate non-reference alleles. These alleles do not have to
    /// be called in any of the samples. Options are base Strings made up of the bases A,C,G,T,N,*, (case insensitive)
    /// or an angle-bracketed ID String ("\<ID\>") or a breakend replacement string as described in the section on
    /// breakends. The '*' allele is reserved to indicate that the allele is missing due to a upstream deletion. If there
    /// are no alternative alleles, then the missing value should be used. Tools processing VCF files are not required
    /// to preserve case in the allele String, except for IDs, which are case sensitive. (String; no whitespace, commas,
    /// or angle-brackets are permitted in the ID String itself)
    pub alt: ArrayRef,

    /// QUAL - quality: Phred-scaled quality score for the assertion made in ALT. i.e. −10log10 prob(call in ALT is
    /// wrong). If ALT is '.' (no variant) then this is −10log10 prob(variant), and if ALT is not '.' this is −10log10
    /// prob(no variant). If unknown, the missing value should be specified. (Numeric)
    pub qual: ArrayRef,

    /// FILTER - filter status: PASS if this position has passed all filters, i.e., a call is made at this position. Otherwise,
    /// if the site has not passed all filters, then a semicolon-separated list of codes for filters that fail. e.g. "q10;s50" might
    /// indicate that at this site the quality is below 10 and the number of samples with data is below 50% of the total
    /// number of samples. '0' is reserved and should not be used as a filter String. If filters have not been applied,
    /// then this field should be set to the missing value. (String, no whitespace or semicolons permitted)
    pub filter: ArrayRef,

    /// INFO - additional information: (String, no whitespace, semicolons, or equals-signs permitted; commas are
    /// permitted only as delimiters for lists of values) INFO fields are encoded as a semicolon-separated series of short
    /// keys with optional values in the format: \<key\>=\<data\>\[,data\]. If no keys are present, the missing value must
    /// be used. Arbitrary keys are permitted, although the following sub-fields are reserved (albeit optional):
    /// - AA: ancestral allele
    /// - AC: allele count in genotypes, for each ALT allele, in the same order as listed
    /// - AF: allele frequency for each ALT allele in the same order as listed: use this when estimated from primary
    ///   data, not called genotypes
    /// - AN: total number of alleles in called genotypes
    /// - BQ: RMS base quality at this position
    /// - CIGAR: cigar string describing how to align an alternate allele to the reference allele
    /// - DB: dbSNP membership
    /// - DP: combined depth across samples, e.g. DP=154
    /// - END: end position of the variant described in this record (for use with symbolic alleles)
    /// - H2: membership in hapmap2
    /// - H3: membership in hapmap3
    /// - MQ: RMS mapping quality, e.g. MQ=52
    /// - MQ0: Number of MAPQ == 0 reads covering this record
    /// - NS: Number of samples with data
    /// - SB: strand bias at this position
    /// - SOMATIC: indicates that the record is a somatic mutation, for cancer genomics
    /// - VALIDATED: validated by follow-up experiment
    /// - 1000G: membership in 1000 Genomes
    pub info: ArrayRef,
    pub samples: HashMap<String, ArrayRef>,
}
