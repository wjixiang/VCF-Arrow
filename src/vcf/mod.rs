//! VCF parsing module.
//!
//! This module provides the core functionality for parsing VCF (Variant Call Format) files
//! and converting them to Apache Arrow format.
//!
//! ## VCF Format Overview
//!
//! VCF is a text format for storing gene sequence variations. A VCF file consists of:
//!
//! 1. **Meta information lines** (starting with `##`): Definitions for contigs, formats, and infos
//! 2. **Header line** (starting with `#`): Column headers including sample names
//! 3. **Data lines**: Each line represents a single variant/call
//!
//! ## Standard VCF Columns
//!
//! | Column | Description |
//! |--------|-------------|
//! | CHROM  | Chromosome identifier |
//! | POS    | Position (1-based) |
//! | ID     | Variant identifier(s) |
//! | REF    | Reference allele |
//! | ALT    | Alternate allele(s) |
//! | QUAL   | Quality score |
//! | FILTER | Filter status |
//! | INFO   | Additional information |
//!
//! ## Specification Reference
//!
//! This module follows the [Variant Call Format (VCF) Version 4.2 Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
//!
//! ## Example
//!
//! ```rust
//! use vcf_arrow::VcfReader;
//!
//! let reader = VcfReader::convert_from_gz("test.vcf.gz")?;
//! let result = reader.parse_into_arrow()?;
//!
//! println!("Chromosome: {:?}", result.chrom);
//! println!("Position: {:?}", result.pos);
//! println!("Samples: {:?}", result.meta.samples);
//! ```

pub mod builders;
pub mod reader;
pub mod types;

pub use builders::*;
pub use reader::*;
pub use types::*;
