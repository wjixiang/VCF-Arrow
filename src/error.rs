//! Error types for VCF parsing operations.
//!
//! This module defines the error types used throughout the VCF-Arrow library.
//! All errors are wrapped in a unified `VcfError` enum for convenient error handling.

/// Errors that can occur during VCF file processing.
///
/// # Variants
///
/// - `LoadVcfError`: Errors encountered while loading or reading VCF files,
///   including file not found, permission denied, or I/O errors.
/// - `ParseVcfError`: Errors encountered while parsing VCF content,
///   including malformed headers, invalid data formats, or type conversion failures.
///
/// # Example
///
/// ```rust
/// use vcf_arrow::{VcfReader, VcfError};
///
/// match VcfReader::convert_from_gz("data.vcf.gz") {
///     Ok(reader) => {
///         match reader.parse_into_arrow() {
///             Ok(result) => { /* use result */ }
///             Err(e) => eprintln!("Parse error: {}", e),
///         }
///     }
///     Err(e) => eprintln!("Load error: {}", e),
/// }
/// ```
#[derive(Debug, thiserror::Error)]
pub enum VcfError {
    /// Error loading VCF file (I/O operations)
    #[error("Load VCF error: {0}")]
    LoadVcfError(String),

    /// Error parsing VCF content (format violations)
    #[error("Parse VCF error: {0}")]
    ParseVcfError(String),
}

impl From<std::io::Error> for VcfError {
    fn from(value: std::io::Error) -> Self {
        VcfError::LoadVcfError(value.to_string())
    }
}
