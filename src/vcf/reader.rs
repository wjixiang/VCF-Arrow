//! VCF reader module.
//!
//! This module provides the `VcfReader` struct, which is the main entry point
//! for loading and parsing VCF files.
//!
//! ## Loading VCF Files
//!
//! `VcfReader` supports three input methods:
//!
//! - **File path**: Load directly from a gzipped `.vcf.gz` file via `convert_from_gz()`
//! - **Bytes**: Load from a `Bytes` buffer containing gzipped VCF data via `convert_from_gz_bytes()`
//! - **String**: Load from an in-memory string via `convert_from_str()`
//!
//! ## Parsing to Arrow
//!
//! After creating a `VcfReader`, use `parse_into_arrow()` to convert the VCF content
//! into Apache Arrow arrays. This produces a `VcfParseResult` containing:
//!
//! - Standard VCF columns as Arrow arrays (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
//! - Metadata (contigs, format definitions, sample names)
//! - Sample data as a map of format IDs to Arrow arrays
//!
//! ## Specification Reference
//!
//! This module follows the [Variant Call Format (VCF) Version 4.2 Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

use std::fs::File;
use std::io::{BufReader, Cursor, Read};
use std::sync::Arc;

use arrow::array::{Int64Builder, StringBuilder};
use bytes::Bytes;
use flate2::read::GzDecoder;

use crate::error::VcfError;
use super::builders::{build_samples, VcfSampleFrameBuilder};
use super::types::{Contig, FormatDef, VcfMeta, VcfParseResult};

/// Core component for converting VCF data into Apache Arrow format.
///
/// `VcfReader` provides a simple interface to load VCF files (either gzipped or plain text)
/// and parse them into Arrow arrays for high-performance data analysis.
///
/// # Example
///
/// ```rust
/// use vcf_arrow::VcfReader;
///
/// // Load from gzipped file
/// let reader = VcfReader::convert_from_gz("variants.vcf.gz")?;
/// let result = reader.parse_into_arrow()?;
///
/// // Access standard columns
/// println!("Chromosomes: {:?}", result.chrom);
/// println!("Positions: {:?}", result.pos);
/// ```
pub struct VcfReader {
    str_content: String,
}

impl VcfReader {
    /// Creates a new `VcfReader` from a gzipped VCF file at the given path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the gzipped VCF file (`.vcf.gz`)
    ///
    /// # Returns
    ///
    /// Returns `Ok(VcfReader)` if the file was successfully read and decompressed,
    /// or `Err(VcfError)` if an I/O error occurred.
    ///
    /// # Example
    ///
    /// ```rust
    /// use vcf_arrow::VcfReader;
    ///
    /// let reader = VcfReader::convert_from_gz("data/test.vcf.gz")?;
    /// ```
    pub fn convert_from_gz(path: &str) -> Result<Self, VcfError> {
        let f = File::open(path)?;
        let buffer = BufReader::new(f);
        let mut gz = GzDecoder::new(buffer);
        let mut content = String::new();
        gz.read_to_string(&mut content)?;

        Ok(Self {
            str_content: content,
        })
    }

    /// Creates a new `VcfReader` from a gzipped byte buffer.
    ///
    /// This is useful when VCF data is already in memory, such as from a network
    /// download or embedded resource.
    ///
    /// # Arguments
    ///
    /// * `bytes` - Gzipped VCF data as `Bytes`
    ///
    /// # Returns
    ///
    /// Returns `Ok(VcfReader)` if the data was successfully decompressed,
    /// or `Err(VcfError)` if an I/O error occurred.
    ///
    /// # Example
    ///
    /// ```rust
    /// use vcf_arrow::VcfReader;
    /// use bytes::Bytes;
    ///
    /// let data = Bytes::from_static(include_bytes!("../tests/test.vcf.gz"));
    /// let reader = VcfReader::convert_from_gz_bytes(data)?;
    /// ```
    pub fn convert_from_gz_bytes(bytes: Bytes) -> Result<Self, VcfError> {
        let cursor = Cursor::new(bytes);
        let mut gz = GzDecoder::new(cursor);
        let mut content = String::new();
        gz.read_to_string(&mut content)?;

        Ok(Self {
            str_content: content,
        })
    }

    /// Creates a new `VcfReader` from a plain string containing VCF content.
    ///
    /// This is useful when the VCF data is already decompressed in memory.
    ///
    /// # Arguments
    ///
    /// * `content` - VCF content as a string slice
    ///
    /// # Returns
    ///
    /// Always returns `Ok(VcfReader)` since no decompression is needed.
    ///
    /// # Example
    ///
    /// ```rust
    /// use vcf_arrow::VcfReader;
    ///
    /// let content = std::fs::read_to_string("data.vcf")?;
    /// let reader = VcfReader::convert_from_str(&content)?;
    /// ```
    pub fn convert_from_str(content: &str) -> Result<Self, VcfError> {
        Ok(Self {
            str_content: content.to_string(),
        })
    }

    /// Parses the VCF content into Apache Arrow arrays.
    ///
    /// This is the main conversion method that transforms VCF data into a
    /// structured `VcfParseResult` containing:
    ///
    /// - **Standard columns**: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    /// - **Metadata**: Contigs, format definitions, and sample names
    /// - **Sample data**: Dynamic arrays keyed by format field IDs
    ///
    /// # Returns
    ///
    /// Returns `Ok(VcfParseResult)` containing all parsed Arrow arrays,
    /// or `Err(VcfError)` if parsing failed.
    ///
    /// # Example
    ///
    /// ```rust
    /// use vcf_arrow::VcfReader;
    ///
    /// let reader = VcfReader::convert_from_gz("test.vcf.gz")?;
    /// let result = reader.parse_into_arrow()?;
    ///
    /// // Access standard columns
    /// let chrom_array = result.chrom;
    /// let pos_array = result.pos;
    ///
    /// // Access sample-specific data
    /// if let Some(gt_array) = result.samples.get("GT") {
    ///     println!("Genotype data: {:?}", gt_array);
    /// }
    /// ```
    pub fn parse_into_arrow(&self) -> Result<VcfParseResult, VcfError> {
        let mut meta = VcfMeta::default();

        let mut chrom_array_builder = StringBuilder::new();
        let mut pos_array_builder = Int64Builder::new();
        let mut id_array_builder = StringBuilder::new();
        let mut ref_array_builder = StringBuilder::new();
        let mut alt_array_builder = StringBuilder::new();
        let mut qual_array_builder = StringBuilder::new();
        let mut filter_array_builder = StringBuilder::new();
        let mut info_array_builder = StringBuilder::new();

        let mut sample_builder = VcfSampleFrameBuilder::default();

        for line in self.str_content.lines() {
            if line.starts_with("##") {
                if let Some(c) = Contig::parse(line) {
                    meta.contigs.push(c);
                }
                if let Some(f) = FormatDef::parse(line) {
                    meta.formats.push(f.clone());
                    sample_builder.add_sample_builder(f)?;
                }
                if let Some(_info) = line.strip_prefix("##INFO=") {
                }
            } else if line.starts_with("#") {
                let cols: Vec<&str> = line.split('\t').collect();
                for col in cols.iter().skip(9) {
                    meta.samples.push(col.to_string());
                }
            } else if !line.is_empty() {
                let row: Vec<&str> = line.split('\t').collect();

                for (index, col) in row.iter().enumerate() {
                    match index {
                        0 => chrom_array_builder.append_value(col),
                        1 => pos_array_builder.append_value(col.parse().unwrap_or(0)),
                        2 => id_array_builder.append_value(col),
                        3 => ref_array_builder.append_value(col),
                        4 => alt_array_builder.append_value(col),
                        5 => qual_array_builder.append_value(col),
                        6 => filter_array_builder.append_value(col),
                        7 => info_array_builder.append_value(col),
                        8 => {
                            let format_data: Vec<&str> = row
                                .get(8)
                                .ok_or(VcfError::ParseVcfError(format!(
                                    "Cannot get format data in '{}'",
                                    line
                                )))?
                                .split(":")
                                .collect();
                            let sample_data: Vec<&str> = row
                                .get(9)
                                .ok_or(VcfError::ParseVcfError(format!(
                                    "Cannot get sample data in '{}'",
                                    line
                                )))?
                                .split(":")
                                .collect();

                            if format_data.len() != sample_data.len() {
                                Err(VcfError::ParseVcfError(format!(
                                    "Formant data length mismatch: {}",
                                    line
                                )))?
                            }

                            for (i, r) in sample_data.iter().enumerate() {
                                let format_id =
                                    format_data.get(i).ok_or(VcfError::ParseVcfError(format!(
                                        "Missing format id at index '{}'",
                                        i
                                    )))?;

                                let builder = sample_builder
                                    .builder_map
                                    .get_mut(*format_id)
                                    .ok_or(VcfError::ParseVcfError(format!(
                                        "No ArrayBuilder registered for id '{}'",
                                        format_id
                                    )))?;

                                builder.arrow_builder.append_from_str(r)?;
                            }
                        }
                        _ => break,
                    }
                }
            }
        }

        let sample_arrow_list = build_samples(&mut sample_builder.builder_map)?;

        Ok(VcfParseResult {
            meta,
            samples: sample_arrow_list,
            chrom: Arc::new(chrom_array_builder.finish()),
            pos: Arc::new(pos_array_builder.finish()),
            id: Arc::new(id_array_builder.finish()),
            _ref: Arc::new(ref_array_builder.finish()),
            alt: Arc::new(alt_array_builder.finish()),
            qual: Arc::new(qual_array_builder.finish()),
            filter: Arc::new(filter_array_builder.finish()),
            info: Arc::new(info_array_builder.finish()),
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    static PATH: &str = "./tests/test.vcf.gz";

    fn get_reader() -> VcfReader {
        VcfReader::convert_from_gz(PATH).unwrap()
    }

    #[test]
    fn test_read_vcf_str() {
        let reader = get_reader();
        let line_vec: Vec<&str> = reader.str_content.lines().collect();
        println!("{:#?}", &line_vec[..10]);
    }

    #[test]
    fn test_read_vcf_meta() {
        let reader = get_reader();
        let res = reader.parse_into_arrow().unwrap();
        assert!(res._ref.as_ref().len() > 0);
        dbg!(res);
    }
}
