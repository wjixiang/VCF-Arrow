//! VCF Arrow builders.
//!
//! This module provides builders for converting VCF data to Apache Arrow arrays.
//! It defines the `DynamicBuilder` trait and implementations for various Arrow array types,
//! as well as the `VcfSampleBuilder` and `VcfSampleFrameBuilder` for handling sample data.
//!
//! ## Dynamic Building
//!
//! The `DynamicBuilder` trait allows runtime selection of appropriate array builders
//! based on VCF format type definitions. This enables automatic type conversion
//! from VCF string data to proper Arrow array types.
//!
//! ## Supported Types
//!
//! | VCF Type | Arrow Builder |
//! |----------|---------------|
//! | String   | StringBuilder |
//! | Integer  | Int32Builder  |
//! | Float    | Float32Builder |

use std::collections::HashMap;

use arrow::array::{ArrayRef, Float32Builder, Int32Builder, StringBuilder};

use super::types::{FormatDef, VcfSample};
use crate::error::VcfError;

/// A trait for building Arrow arrays from string values dynamically.
///
/// This trait enables runtime type conversion from VCF string data to Arrow arrays.
/// Implementations handle parsing string values into the appropriate Arrow array type.
pub trait DynamicBuilder {
    /// Appends a value from a string representation.
    ///
    /// # Arguments
    ///
    /// * `value` - String representation of the value to append
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or `Err(VcfError)` if the value could not be parsed.
    fn append_from_str(&mut self, value: &str) -> Result<(), VcfError>;

    /// Appends a null value to the array.
    fn append_null(&mut self);

    /// Finalizes the builder and returns the constructed Arrow array.
    ///
    /// # Returns
    ///
    /// Returns `Ok(ArrayRef)` containing the finished Arrow array.
    fn build(&mut self) -> Result<ArrayRef, VcfError>;
}

impl DynamicBuilder for StringBuilder {
    fn append_from_str(&mut self, value: &str) -> Result<(), VcfError> {
        self.append_value(value);
        Ok(())
    }

    fn append_null(&mut self) {
        self.append_null();
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

    fn append_null(&mut self) {
        self.append_null();
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

    fn append_null(&mut self) {
        self.append_null();
    }

    fn build(&mut self) -> Result<ArrayRef, VcfError> {
        Ok(std::sync::Arc::new(self.finish()))
    }
}

/// Builder for a single FORMAT field's sample data.
///
/// `VcfSampleBuilder` wraps an Arrow array builder and its corresponding
/// format definition. It handles type-specific value appending for FORMAT fields
/// like GT (genotype), DP (depth), GQ (quality), etc.
pub struct VcfSampleBuilder {
    /// The underlying Arrow array builder
    pub arrow_builder: Box<dyn DynamicBuilder>,

    /// The format definition this builder corresponds to
    format_def: FormatDef,
}

impl VcfSampleBuilder {
    /// Creates a new `VcfSampleBuilder` from a format definition.
    ///
    /// The appropriate Arrow builder is selected based on the format type:
    /// - "String" -> StringBuilder
    /// - "Integer" -> Int32Builder
    /// - "Float" -> Float32Builder
    ///
    /// # Arguments
    ///
    /// * `format_def` - The format definition specifying the field ID, type, and description
    ///
    /// # Returns
    ///
    /// Returns `Ok(VcfSampleBuilder)` if the format type is supported,
    /// or `Err(VcfError)` if the type is not recognized.
    ///
    /// # Example
    ///
    /// ```rust
    /// use vcf_arrow::vcf::types::FormatDef;
    /// use vcf_arrow::vcf::builders::VcfSampleBuilder;
    ///
    /// let format_def = FormatDef {
    ///     id: "GT".to_string(),
    ///     number: "1".to_string(),
    ///     type_: "String".to_string(),
    ///     description: "Genotype".to_string(),
    /// };
    /// let builder = VcfSampleBuilder::new(format_def)?;
    /// ```
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

/// Builder for managing all sample FORMAT data in a VCF file.
///
/// `VcfSampleFrameBuilder` maintains a collection of `VcfSampleBuilder` instances,
/// one for each unique FORMAT field encountered in the VCF header. It provides
/// methods to add new format builders and finalize all sample data into Arrow arrays.
#[derive(Default)]
pub struct VcfSampleFrameBuilder {
    /// Map of FORMAT ID to its corresponding builder
    pub builder_map: HashMap<String, VcfSampleBuilder>,
}

impl VcfSampleFrameBuilder {
    /// Creates a new `VcfSampleFrameBuilder` from a list of format definitions.
    ///
    /// # Arguments
    ///
    /// * `def_vec` - List of format definitions from the VCF header
    ///
    /// # Returns
    ///
    /// Returns `Ok(VcfSampleFrameBuilder)` with builders initialized for each format.
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

    /// Adds a new format builder based on a format definition.
    ///
    /// If a builder for this format ID already exists, it will be replaced.
    ///
    /// # Arguments
    ///
    /// * `format_def` - The format definition to add
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or `Err(VcfError)` if the format type is unsupported.
    pub fn add_sample_builder(&mut self, format_def: FormatDef) -> Result<(), VcfError> {
        self.builder_map
            .insert(format_def.id.clone(), VcfSampleBuilder::new(format_def)?);
        Ok(())
    }
}

/// Finalizes all sample builders into a map of Arrow arrays.
///
/// This function takes the builder map and converts each builder into its
/// final Arrow array, returning a map from FORMAT ID to the resulting array.
///
/// # Arguments
///
/// * `builder_map` - Map of FORMAT ID to VcfSampleBuilder
///
/// # Returns
///
/// Returns `Ok(HashMap<String, ArrayRef>)` containing all sample data arrays,
/// or `Err(VcfError)` if any builder failed to finalize.
///
/// # Example
///
/// ```rust
/// use vcf_arrow::vcf::builders::{VcfSampleFrameBuilder, build_samples};
///
/// let mut frame_builder = VcfSampleFrameBuilder::default();
/// // ... add sample builders ...
/// let sample_arrays = build_samples(&mut frame_builder.builder_map)?;
/// ```
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
