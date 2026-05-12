//! # VCF-Arrow
//!
//! A high-performance VCF (Variant Call Format) parser built with Rust and Apache Arrow.
//!
//! ## Overview
//!
//! VCF-Arrow converts VCF data into Apache Arrow format, enabling seamless integration
//! with Arrow-compatible analysis frameworks such as Polars, DataFusion, and PyArrow.
//!
//! ## Features
//!
//! - **High-performance parsing**: Built with Rust for memory safety and speed
//! - **Arrow integration**: Outputs standard Apache Arrow arrays for easy data manipulation
//! - **VCF 4.2 compliant**: Follows the [VCF Version 4.2 Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
//! - **Gzip support**: Native support for compressed `.vcf.gz` files
//!
//! ## Usage
//!
//! ### Reading a gzipped VCF file
//!
//! ```rust
//! use vcf_arrow::VcfReader;
//!
//! let reader = VcfReader::convert_from_gz("data.vcf.gz")?;
//! let result = reader.parse_into_arrow()?;
//! ```
//!
//! ### Reading from bytes
//!
//! ```rust
//! use vcf_arrow::VcfReader;
//! use bytes::Bytes;
//!
//! let bytes = Bytes::from_static(b"...");
//! let reader = VcfReader::convert_from_gz_bytes(bytes)?;
//! let result = reader.parse_into_arrow()?;
//! ```
//!
//! ### Reading from string
//!
//! ```rust
//! use vcf_arrow::VcfReader;
//!
//! let content = std::fs::read_to_string("data.vcf")?;
//! let reader = VcfReader::convert_from_str(&content)?;
//! let result = reader.parse_into_arrow()?;
//! ```
//!
//! ## VCF Data Structure
//!
//! The conversion result contains:
//! - **Meta information**: Contigs, format definitions, info definitions, and sample names
//! - **Standard columns**: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
//! - **Sample data**: Dynamic arrays keyed by format ID (e.g., GT, DP, GQ)

pub mod error;
pub mod vcf;
