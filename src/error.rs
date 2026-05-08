#[derive(Debug, thiserror::Error)]
pub enum VcfError {
    #[error("Load VCF error: {0}")]
    LoadVcfError(String),

    #[error("Parse VCF error: {0}")]
    ParseVcfError(String),
}

impl From<std::io::Error> for VcfError {
    fn from(value: std::io::Error) -> Self {
        VcfError::LoadVcfError(value.to_string())
    }
}
