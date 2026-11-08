# splice-vus-deconv

**Cell-type-resolved prioritization of non-coding splice and cis-regulatory variants coupling 1Mb AlphaGenome Sashimi junction tracks, simplex-projected NNLS cellular deconvolution, and AlphaFold isoform PAE domain packing.**

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0008--4014--5753-a6ce39.svg)](https://orcid.org/0009-0008-4014-5753)

Developed by the **Nguyen Laboratory for Molecular Genetics, Chromatin Regulation & Oncology** at **Schema Labs Open Science Institute** ([schemalabs.co.uk/research](https://schemalabs.co.uk/research)).

## Scientific Motivation

Over 98% of human ClinVar Variants of Uncertain Significance (VUS) reside in intronic or cis-regulatory regions where bulk tissue RNA-seq masks cell-lineage-specific exon skipping and cryptic splice activation. `splice-vus-deconv` bridges three orthogonal modalities in a lightweight, reproducible Python toolkit:

1. **1Mb Sashimi Junction Perturbation (`splice_vus_deconv.sashimi`)** — Aggregates distance-attenuated exon-inclusion $\Delta\Psi$ across 1Mb cis-regulatory intervals alongside 9-mer U2/U12 5' donor motif log-odds disruption.
2. **Simplex-Projected NNLS Cellular Deconvolution (`splice_vus_deconv.deconv`)** — Resolves cell-type proportions $w \in \Delta^{K-1}$ from bulk transcriptomic mixtures following our benchmark methodology in *Nucleic Acids Research* (2024; `10.1093/nar/gkae267`) and perturbation subtyping in `PINSPlus` (*Bioinformatics* 2019; `10.1093/bioinformatics/bty1049`).
3. **AlphaFold Isoform PAE Domain Destabilization (`splice_vus_deconv.pae`)** — Evaluates inter-domain Pairwise Aligned Error ($\text{PAE}$, Å) shifts across retained residue blocks to distinguish benign alternative isoforms from fold-disrupting cryptic splice events.

## Quickstart

```bash
pip install -e .
python -m unittest discover -s tests -v
splice-vus-deconv --benchmark
```

## Citation

- Nguyen, H., Nguyen, H., Tran, D., Draghici, S., & Nguyen, T. (2024). *Fourteen years of cellular deconvolution: methodology, applications, technical evaluation and outstanding challenges.* **Nucleic Acids Research**, 52(9), 4761–4783. https://doi.org/10.1093/nar/gkae267
- Nguyen, H., Shrestha, S., Draghici, S., & Nguyen, T. (2019). *PINSPlus: a tool for tumor subtype discovery in integrated genomic data.* **Bioinformatics**, 35(16), 2843–2846. https://doi.org/10.1093/bioinformatics/bty1049
