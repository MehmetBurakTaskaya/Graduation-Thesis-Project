# Upstream sequence analysis of bacterial housekeeping genes

Analysis code for my B.Sc. graduation project in Molecular Biology and Genetics at Acıbadem University (2023), supervised by Dr. Erkan Mozioğlu.

**Thesis title:** *Decoding Sequence Patterns and Conserved Motifs in Housekeeping Gene Upstream Regions: Unveiling Deeper Implications in Escherichia coli*

## What the code does

The script compares the 33-nucleotide regions directly upstream of bacterial genes and looks for shared sequence patterns. It:

- checks which upstream sequences are identical;
- builds UPGMA trees from identity distances (Biopython);
- counts nucleotide frequencies at each position, as line and stacked-bar plots;
- finds the most frequent subsequence of every length, across all sequences or per gene;
- searches each sequence for ribosome-binding-site (Shine–Dalgarno-like) motifs and records the motif found and its position. Two motif sets are used: a canonical set and a selected set.

## Data (`Datasets/`)

| File | Contents |
|---|---|
| `sequences.xlsx` | 300 upstream sequences of six housekeeping genes (*rpoB*, *recA*, *infB*, *gyrB*, *gapA*, *dnaE*), 50 per gene, from 222 bacterial species and strains |
| `Escherichia_coli_sequences.xlsx` | 49 upstream sequences of 48 genes from *Escherichia coli* str. K-12 substr. MG1655 |

Both files have the columns `ID`, `GeneName`, `SpecName`, `Pre-Sequence` (the upstream sequence) and `Pre-Length` (33 for all rows).

## Running it

The code is one script, `GraduationThesisProject.py`, split into `# region` blocks. Each block is one analysis and is meant to be run on its own, for example cell by cell in VS Code or Spyder.

Requirements: Python 3 with `pandas`, `openpyxl`, `numpy`, `matplotlib` and `biopython`.

```bash
pip install pandas openpyxl numpy matplotlib biopython
```

Before running:

- **Set the paths.** The script sets its working directory with hard-coded Windows paths at the top and in the "Gene Phylogenetic Tree" block. Change them to your own folder, or run from inside `Datasets/`.
- **Skip or adapt the "Fasta Downloader" block.** It uses a placeholder URL (`example.com`) and does not download anything as written.
- **Some blocks write files.** `RBS_results.xlsx` and the `Fasta_files/` or `Fasta_sequences/` folders are written to the working directory.

## Summary (adapted from the thesis abstract)

This graduation thesis project explores sequence patterns and conserved motifs in the upstream regions of housekeeping genes. Several housekeeping gene sequences are compared within a dataset to identify conserved areas and common motifs. The nucleotide sequences of the housekeeping genes rpoB, recA, infB, gyrB, gapA and dnaE are examined both individually and as a combined dataset. The findings are then applied to gene sequences from *Escherichia coli* str. K-12 substr. MG1655 and compared with earlier outcomes.

## Author

Mehmet Burak Taşkaya · [ORCID 0009-0002-0178-5199](https://orcid.org/0009-0002-0178-5199)
