# VCF-Arrow

A high-performance VCF (Variant Call Format) parser built with Rust and Apache Arrow, designed for modern bioinformatics data analysis pipelines.

## Features

- **Complete Raw Data Parsing** — Fully parses VCF metadata (contigs, FORMAT, INFO) and data body with zero information loss
- **Automatic Schema Conversion** — Automatically generates typed Arrow arrays for sample data based on FORMAT definitions, producing highly structured output that reduces downstream data cleaning overhead
- **Arrow Ecosystem Integration** — All data is converted to Arrow `ArrayRef`, enabling seamless integration with Arrow-compatible analysis frameworks (e.g., Polars, DataFusion, PyArrow) without additional serialization/deserialization
- **gzip Support** — Native `.vcf.gz` decompression via `flate2`

## Supported Types

### Data Types

| VCF Type  | Arrow Type     |
| --------- | -------------- |
| `Integer` | `Int32Array`   |
| `Float`   | `Float32Array` |
| `String`  | `StringArray`  |

### Meta-information

- [x] Contig
- [x] Filter field format
- [x] Sample field format
- [] File format
- [] Information field format
- [] Alternative allele field format
- [] Assembly field format
- [] Pedigree field format

## Usage

Add `vcf-arrow` to your `Cargo.toml`:

```toml
[dependencies]
vcf-arrow = "0.1.2"
```

### Parse a `.vcf.gz` file

```rust
use vcf_arrow::vcf::VcfReader;

let reader = VcfReader::load_gz("sample.vcf.gz")?;
let result = reader.parse_into_arrow()?;

println!("Samples: {:?}", result.meta.samples);
println!("CHROM:  {:?}", result.chrom);
println!("POS:    {:?}", result.pos);
println!("Sample data: {:?}", result.samples);
```

### Parse from a string

```rust
use vcf_arrow::vcf::VcfReader;

let reader = VcfReader::convert_from_str(include_str!("sample.vcf"))?;
let result = reader.parse_into_arrow()?;
```

## Core Structs

| Struct           | Description                                                                                                                                              |
| ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `VcfReader`      | Entry point for loading and parsing VCF files                                                                                                            |
| `VcfParseResult` | Parsed result containing metadata (`VcfMeta`), fixed columns (`chrom`, `pos`, `id`, `alt`, `qual`, `filter`, `info`), and sample data (`Vec<VcfSample>`) |
| `VcfMeta`        | Parsed metadata: contigs, FORMAT definitions, INFO definitions, and sample names                                                                         |
| `VcfSample`      | A single sample field with its `FormatDef` and the corresponding Arrow `ArrayRef`                                                                        |

## Reference

This library follows the [VCF Version 4.2 Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

## License

Licensed under either of [MIT](https://opensource.org/licenses/MIT) or [Apache-2.0](https://www.apache.org/licenses/LICENSE-2.0) at your option.
