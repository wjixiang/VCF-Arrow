//! VCF reader module
//! This module use [ Variant Call Format(VCF) Version 4.2 Specification ](https://samtools.github.io/hts-specs/VCFv4.2.pdf) as standard reference

use std::fs::File;
use std::io::{BufReader, Cursor, Read};
use std::sync::Arc;

use arrow::array::{Int64Builder, StringBuilder};
use bytes::Bytes;
use flate2::read::GzDecoder;

use crate::error::VcfError;
use super::builders::{build_samples, VcfSampleFrameBuilder};
use super::types::{Contig, FormatDef, VcfMeta, VcfParseResult};

/// Core component that capable to convert VCF data into Appache-Arrow
pub struct VcfReader {
    str_content: String,
}

impl VcfReader {
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

    pub fn convert_from_gz_bytes(bytes: Bytes) -> Result<Self, VcfError> {
        let cursor = Cursor::new(bytes);
        let mut gz = GzDecoder::new(cursor);
        let mut content = String::new();
        gz.read_to_string(&mut content)?;

        Ok(Self {
            str_content: content,
        })
    }

    /// Core conversion method
    pub fn convert_from_str(content: &str) -> Result<Self, VcfError> {
        Ok(Self {
            str_content: content.to_string(),
        })
    }

    pub fn parse_into_arrow(&self) -> Result<VcfParseResult, VcfError> {
        let mut meta = VcfMeta::default();

        // Builders for each column
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
                // INFO fields are parsed but not used for now
                if let Some(_info) = line.strip_prefix("##INFO=") {
                    // meta.infos.push(info);
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
                        0 => chrom_array_builder.append_value(col), // #CHROM
                        1 => pos_array_builder.append_value(col.parse().unwrap_or(0)), // POS
                        2 => id_array_builder.append_value(col),    // ID
                        3 => ref_array_builder.append_value(col),   // REF
                        4 => alt_array_builder.append_value(col),   // ALT
                        5 => qual_array_builder.append_value(col),  // QUAL
                        6 => filter_array_builder.append_value(col), // FILTER
                        7 => info_array_builder.append_value(col),  // INFO
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
                                // Check the length of both format column and sample column
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
                        _ => break, // Skip sample columns
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
