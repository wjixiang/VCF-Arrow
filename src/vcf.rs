//! VCF parsing module
//! This module use [ Variant Vall Format(VCF) Version 4.2 Specification ](https://samtools.github.io/hts-specs/VCFv4.2.pdf) as standard reference

use std::collections::HashMap;
use std::io::{Cursor, Read};
use std::sync::Arc;
use std::{fs::File, io::BufReader};

use arrow::array::{ArrayRef, Float32Builder, Int32Builder, Int64Builder, StringBuilder};
use bytes::Bytes;
use flate2::read::GzDecoder;

use crate::error::VcfError;

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
    fn parse(line: &str) -> Option<Self> {
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
    fn parse(line: &str) -> Option<Self> {
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

/// Core component that capable to convert VCF data into Appache-Arrow
pub struct VcfReader {
    str_content: String,
}

#[derive(Debug)]
pub struct VcfSample {
    pub format_def: FormatDef,
    pub array: ArrayRef,
}

pub trait DynamicBuilder {
    fn append_from_str(&mut self, value: &str) -> Result<(), VcfError>;
    fn build(&mut self) -> Result<ArrayRef, VcfError>;
}

impl DynamicBuilder for StringBuilder {
    fn append_from_str(&mut self, value: &str) -> Result<(), VcfError> {
        self.append_value(value);

        Ok(())
    }

    fn build(&mut self) -> Result<ArrayRef, VcfError> {
        Ok(Arc::new(self.finish()))
    }
}

impl DynamicBuilder for Float32Builder {
    fn append_from_str(&mut self, value: &str) -> Result<(), VcfError> {
        self.append_value(match value.parse() {
            Ok(it) => it,
            Err(err) => {
                return Err(VcfError::ParseVcfError(format!(
                    "Cannot convert '{}' into f32, {}",
                    value, err
                )));
            }
        });
        Ok(())
    }

    fn build(&mut self) -> Result<ArrayRef, VcfError> {
        Ok(Arc::new(self.finish()))
    }
}

impl DynamicBuilder for Int32Builder {
    fn append_from_str(&mut self, value: &str) -> Result<(), VcfError> {
        self.append_value(value.parse().map_err(|err| {
            VcfError::ParseVcfError(format!("Cannot convert '{}' into int32, {}", value, err))
        })?);
        Ok(())
    }

    fn build(&mut self) -> Result<ArrayRef, VcfError> {
        Ok(Arc::new(self.finish()))
    }
}

pub struct VcfSampleBuilder {
    pub arrow_builder: Box<dyn DynamicBuilder>,
    format_def: FormatDef,
}

impl VcfSampleBuilder {
    pub fn new(format_def: FormatDef) -> Result<Self, VcfError> {
        let array_builder: Box<dyn DynamicBuilder> = match format_def.type_.as_str() {
            "Float" => Box::new(Float32Builder::new()),
            "Integer" => Box::new(Int32Builder::new()),
            "String" => Box::new(StringBuilder::new()),
            _ => Err(VcfError::ParseVcfError(format!(
                "Unsupported vcf field type: {}",
                format_def.type_
            )))?,
        };
        Ok(Self {
            arrow_builder: array_builder,
            format_def,
        })
    }
}

#[derive(Default)]
pub struct VcfSampleFrameBuilder {
    pub builder_map: HashMap<String, VcfSampleBuilder>,
}

impl VcfSampleFrameBuilder {
    pub fn new(def_vec: Vec<FormatDef>) -> Result<Self, VcfError> {
        let mut smaple_container: HashMap<String, VcfSampleBuilder> = HashMap::new();

        for i in def_vec {
            let id = i.id.clone();
            let builder = VcfSampleBuilder::new(i)?;
            smaple_container.insert(id, builder);
        }

        Ok(Self {
            builder_map: smaple_container,
        })
    }

    pub fn add_sample_builder(&mut self, format_def: FormatDef) -> Result<(), VcfError> {
        self.builder_map
            .insert(format_def.id.clone(), VcfSampleBuilder::new(format_def)?);
        Ok(())
    }
}

/// VCF-to-Arrow conversion result struct
pub struct VcfParseResult {
    pub meta: VcfMeta,
    pub chrom: ArrayRef,
    pub pos: ArrayRef,
    pub id: ArrayRef,
    pub alt: ArrayRef,

    /// Phred quality score
    pub qual: ArrayRef,
    pub filter: ArrayRef,
    pub info: ArrayRef,
    pub samples: Vec<VcfSample>,
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

        let sample_arrow_list: Vec<VcfSample> = sample_builder
            .builder_map
            .iter_mut()
            .map(|(_index, builder)| {
                Ok(VcfSample {
                    array: builder.arrow_builder.build()?,
                    format_def: builder.format_def.clone(),
                })
            })
            .collect::<Result<Vec<VcfSample>, VcfError>>()?;

        Ok(VcfParseResult {
            meta,
            samples: sample_arrow_list,
            chrom: Arc::new(chrom_array_builder.finish()),
            pos: Arc::new(pos_array_builder.finish()),
            id: Arc::new(id_array_builder.finish()),
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
        dbg!(res.samples.first().unwrap().array.clone());
    }
}
