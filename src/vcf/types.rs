//! VCF data types.
//!
//! This module defines the core data structures used to represent VCF metadata and
//! parsed results. It includes types for contigs, format definitions, info definitions,
//! and the final VCF-to-Arrow conversion result.
//!
//! ## Specification Reference
//!
//! This module follows the [Variant Call Format (VCF) Version 4.2 Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

use std::collections::HashMap;

use arrow::array::ArrayRef;

/// Represents a contig definition from the VCF header.
///
/// Contigs define the reference sequences (chromosomes) used in the VCF file.
/// Each contig has an ID, optional length, and optional assembly identifier.
#[derive(Debug, Default)]
pub struct Contig {
    /// Contig identifier (e.g., "chr1", "1", "MT")
    pub id: String,

    /// Length of the contig in bases (optional)
    pub length: Option<u64>,

    /// Assembly identifier (optional)
    pub assembly: Option<String>,
}

impl Contig {
    /// Parses a contig line from the VCF meta section.
    ///
    /// Expects a line in the format: `##contig=<ID=chr1,length=249250621,assembly=GRCh37>`
    ///
    /// # Arguments
    ///
    /// * `line` - A line starting with `##contig=`
    ///
    /// # Returns
    ///
    /// Returns `Some(Contig)` if parsing succeeded, or `None` if the line format was invalid.
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

/// Represents a FORMAT field definition from the VCF header.
///
/// FORMAT fields define the structure of genotype data for each sample.
/// Common FORMAT fields include GT (genotype), DP (depth), GQ (genotype quality), etc.
#[derive(Debug, Clone, Default)]
pub struct FormatDef {
    /// Format field identifier (e.g., "GT", "DP", "GQ")
    pub id: String,

    /// Number of values (can be a number or "G", "R", "A" per VCF spec)
    pub number: String,

    /// Data type (e.g., "String", "Integer", "Float")
    pub type_: String,

    /// Human-readable description of the field
    pub description: String,
}

impl FormatDef {
    /// Parses a FORMAT definition line from the VCF meta section.
    ///
    /// Expects a line in the format: `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`
    ///
    /// # Arguments
    ///
    /// * `line` - A line starting with `##FORMAT=`
    ///
    /// # Returns
    ///
    /// Returns `Some(FormatDef)` if parsing succeeded, or `None` if the line format was invalid.
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

/// Represents an INFO field definition from the VCF header.
///
/// INFO fields provide additional metadata about variants.
/// Common INFO fields include DP (depth), AF (allele frequency), ANN (annotations), etc.
#[derive(Debug, Default)]
pub struct InfoDef {
    /// Info field identifier
    pub id: String,

    /// Number of values
    pub number: String,

    /// Data type
    pub type_: String,

    /// Human-readable description
    pub description: String,
}

/// VCF metadata container.
///
/// Holds all meta information parsed from the VCF header, including
/// contig definitions, FORMAT definitions, INFO definitions, and sample names.
#[derive(Debug, Default)]
pub struct VcfMeta {
    /// Contig definitions from the header
    pub contigs: Vec<Contig>,

    /// FORMAT field definitions
    pub formats: Vec<FormatDef>,

    /// INFO field definitions
    pub infos: Vec<InfoDef>,

    /// Sample names from the header line
    pub samples: Vec<String>,
}

/// Represents a single sample's data array with its format definition.
///
/// This pairs a format definition (e.g., GT, DP) with the corresponding
/// Arrow array containing the sample data for that format field.
#[derive(Debug)]
pub struct VcfSample {
    /// The format definition for this sample data
    pub format_def: FormatDef,

    /// The Arrow array containing sample values
    pub array: ArrayRef,
}

/// VCF-to-Arrow conversion result.
///
/// This is the primary output structure from `VcfReader::parse_into_arrow()`.
/// It contains all standard VCF columns as Arrow arrays, metadata, and
/// sample-specific data.
///
/// # Standard Columns
///
/// The first 8 columns follow the VCF specification:
///
/// | Field | Type | Description |
/// |-------|------|-------------|
/// | chrom | String | Chromosome/contig identifier |
/// | pos | Int64 | 1-based position |
/// | id | String | Variant identifier(s) |
/// | ref | String | Reference allele |
/// | alt | String | Alternate allele(s) |
/// | qual | String | Phred-scaled quality score |
/// | filter | String | Filter status |
/// | info | String | Additional information |
///
/// # Sample Data
///
/// Sample-specific data is stored in the `samples` HashMap, keyed by FORMAT ID.
/// For example, if the VCF contains `##FORMAT=<ID=GT,...>` and `##FORMAT=<ID=DP,...>`,
/// the samples map will contain entries for "GT" and "DP".
#[derive(Debug)]
pub struct VcfParseResult {
    /// VCF metadata (contigs, formats, samples)
    pub meta: VcfMeta,

    /// CHROM - chromosome identifier.
    ///
    /// An identifier from the reference genome or an angle-bracketed ID String ("\<ID\>")
    /// pointing to a contig in the assembly file. All entries for a specific
    /// CHROM should form a contiguous block within the VCF file.
    ///
    /// (String, no whitespace permitted, Required).
    pub chrom: ArrayRef,

    /// POS - reference position.
    ///
    /// The reference position, with the 1st base having position 1. Positions are sorted numerically,
    /// in increasing order, within each reference sequence CHROM. It is permitted to have multiple
    /// records with the same POS. Telomeres are indicated by using positions 0 or N+1, where N is
    /// the length of the corresponding chromosome or contig.
    ///
    /// (Integer, Required)
    pub pos: ArrayRef,

    /// ID - variant identifier.
    ///
    /// Semicolon-separated list of unique identifiers where available. If this is a dbSNP variant
    /// it is encouraged to use the rs number(s). No identifier should be present in more than one
    /// data record. If there is no identifier available, then the missing value should be used.
    ///
    /// (String, no whitespace or semicolons permitted)
    pub id: ArrayRef,

    /// REF - reference allele.
    ///
    /// Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted.
    /// The value in the POS field refers to the position of the first base in the String.
    /// For simple insertions and deletions in which either the REF or one of the ALT alleles
    /// would otherwise be null/empty, the REF and ALT Strings must include the base before the
    /// event (which must be reflected in the POS field).
    ///
    /// (String, Required)
    pub _ref: ArrayRef,

    /// ALT - alternate allele(s).
    ///
    /// Comma separated list of alternate non-reference alleles. These alleles do not have to
    /// be called in any of the samples. Options are base Strings made up of the bases A,C,G,T,N,
    /// *, (case insensitive) or an angle-bracketed ID String ("\<ID\>") or a breakend replacement
    /// string. The '*' allele is reserved to indicate that the allele is missing due to a upstream
    /// deletion. If there are no alternative alleles, then the missing value should be used.
    ///
    /// (String; no whitespace, commas, or angle-brackets are permitted in the ID String itself)
    pub alt: ArrayRef,

    /// QUAL - quality score.
    ///
    /// Phred-scaled quality score for the assertion made in ALT. i.e. −10log10 prob(call in ALT
    /// is wrong). If ALT is '.' (no variant) then this is −10log10 prob(variant), and if ALT is
    /// not '.' this is −10log10 prob(no variant). If unknown, the missing value should be specified.
    ///
    /// (Numeric)
    pub qual: ArrayRef,

    /// FILTER - filter status.
    ///
    /// PASS if this position has passed all filters, i.e., a call is made at this position.
    /// Otherwise, if the site has not passed all filters, then a semicolon-separated list of
    /// codes for filters that fail. e.g. "q10;s50" might indicate that at this site the quality
    /// is below 10 and the number of samples with data is below 50% of the total number of samples.
    /// If filters have not been applied, then this field should be set to the missing value.
    ///
    /// (String, no whitespace or semicolons permitted)
    pub filter: ArrayRef,

    /// INFO - additional information.
    ///
    /// INFO fields are encoded as a semicolon-separated series of short keys with optional values
    /// in the format: \<key\>=\<data\>\[,data\]. If no keys are present, the missing value must
    /// be used. Arbitrary keys are permitted, although the following sub-fields are reserved:
    ///
    /// - AA: ancestral allele
    /// - AC: allele count in genotypes, for each ALT allele
    /// - AF: allele frequency for each ALT allele
    /// - AN: total number of alleles in called genotypes
    /// - DP: combined depth across samples
    /// - MQ: RMS mapping quality
    /// - NS: Number of samples with data
    /// - SOD: strand bias
    ///
    /// (String, no whitespace, semicolons, or equals-signs permitted)
    pub info: ArrayRef,

    /// Sample data as a map of FORMAT ID to Arrow arrays.
    ///
    /// Each key corresponds to a FORMAT field ID (e.g., "GT", "DP", "GQ"),
    /// and the value is an Arrow array containing the data for all samples
    /// for that format field.
    ///
    /// # Example
    ///
    /// ```rust
    /// // Access genotype data
    /// if let Some(gt_array) = result.samples.get("GT") {
    ///     // Process genotype data
    /// }
    ///
    /// // Access depth data
    /// if let Some(dp_array) = result.samples.get("DP") {
    ///     // Process depth data
    /// }
    /// ```
    pub samples: HashMap<String, ArrayRef>,
}
