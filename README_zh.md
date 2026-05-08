# VCF-Arrow

基于Rust与Apache Arrow的VCF解析工具，使用现代化的数据分析生态处理生信分析

- 完整解析原始数据：能够完全解析元数据与数据主体，提高解析速度的同时不遗漏信息
- 自动Schema转换：结合Format数据自动完成Sample数据的Arrow Schema生成，结果相比于原始VCF数据高度结构化，降低数据清洗负担
- 对接Arrow生态：所有数据统一转换为Arrow的ArrayRef，后续分析可灵活对接多种分析框架，无需额外的序列化/反序列化
- 零氛围式编程：解析结果安全可控

A VCF data parser using Appache Arrow as standard format
