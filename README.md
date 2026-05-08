# VCF-Arrow

A VCF parsing tool based on Rust and Apache Arrow, leveraging modern data analysis ecosystems for bioinformatics processing.

- **Complete Raw Data Parsing**: Fully parses metadata and data body, improving parsing speed without losing information
- **Automatic Schema Conversion**: Automatically generates Arrow Schema for sample data combined with Format data, resulting in highly structured data compared to raw VCF, reducing data cleaning overhead
- **Arrow Ecosystem Integration**: All data is uniformly converted to Arrow's ArrayRef, enabling flexible integration with multiple analysis frameworks for downstream analysis without additional serialization/deserialization
