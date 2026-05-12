# VCF-Arrow

基于 Rust 与 Apache Arrow 的高性能 VCF (Variant Call Format) 解析库，专为现代生信数据分析流程设计。

## 特性

- **完整解析原始数据** — 完全解析 VCF 元数据（contig、FORMAT、INFO）与数据主体，零信息丢失
- **自动 Schema 转换** — 根据 FORMAT 定义自动生成带类型的 Arrow 数组，输出高度结构化数据，降低下游数据清洗负担
- **无缝对接 Arrow 生态** — 所有数据统一转换为 Arrow `ArrayRef`，可无缝对接多种 Arrow 兼容分析框架（如 Polars、DataFusion、PyArrow），无需额外序列化/反序列化
- **gzip 支持** — 通过 `flate2` 原生支持 `.vcf.gz` 文件解压

## 支持的类型映射

| VCF 类型  | Arrow 类型     |
|-----------|---------------|
| `Integer` | `Int32Array`   |
| `Float`   | `Float32Array` |
| `String`  | `StringArray`  |

## 使用方法

在 `Cargo.toml` 中添加依赖：

```toml
[dependencies]
vcf-arrow = "0.1.0"
```

### 解析 `.vcf.gz` 文件

```rust
use vcf_arrow::vcf::VcfReader;

let reader = VcfReader::load_gz("sample.vcf.gz")?;
let result = reader.parse_into_arrow()?;

println!("样本: {:?}", result.meta.samples);
println!("CHROM: {:?}", result.chrom);
println!("POS:   {:?}", result.pos);
println!("样本数据: {:?}", result.samples);
```

### 从字符串解析

```rust
use vcf_arrow::vcf::VcfReader;

let reader = VcfReader::convert_from_str(include_str!("sample.vcf"))?;
let result = reader.parse_into_arrow()?;
```

## 核心结构体

| 结构体 | 说明 |
|--------|------|
| `VcfReader` | 加载和解析 VCF 文件的入口 |
| `VcfParseResult` | 解析结果，包含元数据（`VcfMeta`）、固定列（`chrom`、`pos`、`id`、`alt`、`qual`、`filter`、`info`）和样本数据（`Vec<VcfSample>`） |
| `VcfMeta` | 解析后的元数据：contig、FORMAT 定义、INFO 定义和样本名称 |
| `VcfSample` | 单个样本字段，包含其 `FormatDef` 和对应的 Arrow `ArrayRef` |

## 参考规范

本库遵循 [VCF Version 4.2 规范](https://samtools.github.io/hts-specs/VCFv4.2.pdf)。

## 许可证

基于 [MIT](https://opensource.org/licenses/MIT) 或 [Apache-2.0](https://www.apache.org/licenses/LICENSE-2.0) 许可证，可任选其一。
