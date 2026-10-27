---
title: Cross-Species Animal Alignments Dataset (animals.meg)
type: source
category: dataset
tags: [dataset, animal-alignments, mega-format, cross-species, multiple-sequence-alignment, felsenstein-bootstrap]
created: 2026-08-17
updated: 2026-10-28
raw_source: "data/raw/animal_alignments/animals.meg"
---

# Source: Cross-Species Animal Alignments Dataset

## Summary
The `data/raw/animal_alignments/` directory contains cross-species nucleotide alignment benchmarks stored in native MEGA formats (`.meg`, `.mas`, `.mdsx`, `.mtsx`). It provides an 839-site multi-taxa dataset spanning Mammalia (*Homo sapiens*, *Canis lupus familiaris*, *Mus musculus*, *Rattus norvegicus*, *Panthera pardus*), Reptilia (*Chelonia mydas*), Actinopterygii (*Saurida tumbil*), and divergent bacterial/pathogen outgroups (*Jeotgalicoccus halotolerans*, *Snake plant crown rot tissues*).

## Format Specifications
- **MEGA Header**: `#MEGA`, `!Title animals;`, `DataType=Nucleotide`, `NSeqs=10 NSites=839`.
- **Character Conventions**: Identical base `.`, Missing data `?`, Alignment indel/gap `-`.
- **Associated MEGA Files**:
  - `animals.mas`: MEGA Alignment Session metadata.
  - `animals.mdsx`: MEGA Distance Session matrix.
  - `animals.mtsx`: MEGA Tree Session topology data.

## Reconstructed Phylogeny & Benchmark Results
Using [`scripts/phylo_pipeline.py`](file:///Users/macbookairm2/Documents/GitHub/SuleimanMEGAMolecularEvolutionaryTree/scripts/phylo_pipeline.py), the 9 valid taxa were parsed and reconstructed under both K2P and JC69 substitution models with 100 Felsenstein bootstrap replicates:

- **Mammalian Core Group**: *Homo sapiens* and *Canis lupus familiaris* cluster with robust bootstrap support (71.0%), exhibiting lower pairwise divergence ($d_{\text{K2P}} = 1.0978$).
- **Deep Evolutionary Saturation**: Outgroups such as *Jeotgalicoccus halotolerans* and deeply divergent reptile/fish lineages exceed $d > 1.8$, demonstrating classic transition/transversion substitution saturation on long branches.
- **Artifact Outputs**:
  - `results/animals_analysis/animals_nj_tree_k2p.nwk`
  - `results/animals_analysis/animals_nj_tree_k2p.png`
  - `results/animals_analysis/animals_nj_tree_jc69.nwk`
  - `results/animals_analysis/animals_nj_tree_jc69.png`

## Related Graph Nodes
- Entities: [[Homo_sapiens]], [[Mus_musculus]], [[Rattus_norvegicus]], [[MEGA12]]
- Concepts: [[Phylogenetic_Reconstruction]], [[Substitution_Models]], [[Pairwise_Evolutionary_Distance]], [[Bootstrap_Resampling_and_Phylogenetic_Support]]
- Sources: [[Phylo_Pipeline_Script]]
