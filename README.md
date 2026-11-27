# Bioinformatics Portfolio

Small, reproducible bioinformatics projects built with Python, using public data.
Each project has its own folder with code, a README, and figures.

**Bismillah Iqbal** - Biotechnology | Bioinformatics | Python | Genomics

## Projects

| # | Project | What it does | Tools |
|---|---|---|---|
| 01 | [Sequence Analysis Toolkit](01_sequence_toolkit) | Finds open reading frames, GC content and codon usage for any GenBank sequence. Tested on human TP53 (longest ORF = 393 aa, matching the p53 protein). | Python, Biopython, Matplotlib |
| 02 | [RNA-seq Differential Expression: Lung Cancer](02_rnaseq_de) | Paired tumor vs normal analysis of 19 NSCLC patients (GEO: GSE81089). Finds 14,224 significant genes; normal lung genes such as *SFTPC* and *AGER* are lost, while proliferation and invasion genes such as *TOP2A* and *MMP11* are gained. | Python, PyDESeq2, pandas, seaborn |

## Skills shown
- Sequence analysis and ORF finding
- RNA-seq differential expression (paired design, normalization, multiple-testing correction)
- Data visualization: PCA, volcano plot, MA plot, clustered heatmap
- Working with public data (NCBI GenBank, GEO) and reproducible code

## License
MIT
