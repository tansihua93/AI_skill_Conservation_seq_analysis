# Nucleic_acid_Conservation_seq_analysis
End-to-end workflow for batch nucleotide conservation analysis: from keyword-based retrieval on NCBI, through multiple sequence alignment, to sliding-window conservation scoring, conserved-region extraction, and a publication-style multi-panel figure.

# Skill 结构
ncbi-conservation-analysis/
├── SKILL.md                              # 主入口，描述触发词与工作流
├── scripts/
│   ├── run_pipeline.py                   # 主流程编排器
│   ├── download_ncbi.py                  # NCBI Entrez 下载（限速/重试）
│   ├── sliding_window_identity.py        # 逐核苷酸 identity + 50bp 滑动窗口
│   ├── extract_conserved.py              # 合并相邻窗口、导出 FASTA
│   └── plot_conservation.py              # 4 面板 matplotlib 论文级图
├── references/
│   ├── sars-cov-2-annotation.json        # 默认 SARS-CoV-2 ORF 注释
│   ├── figure_layout.md                  # 4 面板布局规范
│   └── dependencies.md                   # 依赖（biopython/matplotlib/MAFFT）
└── assets/
    ├── config_template.json              # 开箱即用配置
    ├── USAGE.md                          # 运行示例
    └── sample_sequences.fasta            # 测试数据
