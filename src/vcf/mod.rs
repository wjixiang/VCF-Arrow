//! VCF parsing module
//! This module use [ Variant Call Format(VCF) Version 4.2 Specification ](https://samtools.github.io/hts-specs/VCFv4.2.pdf) as standard reference

pub mod builders;
pub mod reader;
pub mod types;

pub use builders::*;
pub use reader::*;
pub use types::*;
