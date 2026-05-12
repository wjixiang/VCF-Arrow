//! VCF Arrow builders
//! Builders for converting VCF data to Apache Arrow arrays

use std::collections::HashMap;

use arrow::array::{ArrayRef, Float32Builder, Int32Builder, StringBuilder};

use super::types::{FormatDef, VcfSample};
use crate::error::VcfError;

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
        Ok(std::sync::Arc::new(self.finish()))
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
        Ok(std::sync::Arc::new(self.finish()))
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
        Ok(std::sync::Arc::new(self.finish()))
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

pub fn build_samples(
    builder_map: &mut HashMap<String, VcfSampleBuilder>,
) -> Result<HashMap<String, ArrayRef>, VcfError> {
    let mut res_map: HashMap<String, ArrayRef> = HashMap::new();
    for (index, builder) in builder_map {
        let array = builder.arrow_builder.build()?;
        res_map.insert(index.clone(), array);
    }

    Ok(res_map)
}
