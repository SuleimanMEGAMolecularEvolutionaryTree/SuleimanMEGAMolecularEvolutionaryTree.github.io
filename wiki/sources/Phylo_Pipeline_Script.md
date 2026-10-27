---
title: Pure-Python Molecular Evolutionary Genetics Pipeline (phylo_pipeline.py)
type: source
category: script
tags: [source, python, algorithms, jc69, k2p, neighbor-joining, bootstrap, felsenstein, pure-python, cli]
created: 2026-08-17
updated: 2026-10-28
raw_source: "scripts/phylo_pipeline.py"
---

# Source: Pure-Python Molecular Evolutionary Genetics Pipeline

## Summary
[`scripts/phylo_pipeline.py`](file:///Users/macbookairm2/Documents/GitHub/SuleimanMEGAMolecularEvolutionaryTree/scripts/phylo_pipeline.py) is a comprehensive, standalone Python pipeline for evolutionary distance estimation, Neighbor-Joining tree construction, and Felsenstein non-parametric bootstrap resampling authored by Suleyman Hajizadeh for a Cambridge MPhil portfolio. It implements mathematical definitions of [[Pairwise_Evolutionary_Distance|JC69 and K2P]] models, Saitou & Nei (1987) [[Neighbor_Joining_Algorithm|Neighbor-Joining]], and Felsenstein (1985) [[Bootstrap_Resampling_and_Phylogenetic_Support|Bootstrap Resampling]] from first principles using only standard numerical libraries (`numpy`, `matplotlib`, `biopython`), completely bypassing GUI dependencies like [[MEGA12]].

## Architecture & Multi-Dataset Support

```
1. Multi-Format Ingest ──> 2. Pairwise/MSA Count ──> 3. Analytical Distance
   (NCBI / .meg / .fas)       (Purines/Pyrimidines)     (JC69 / K2P Matrix)
                                                             │
                                                             ▼
6. Dual Output ◄────── 5. Bootstrap Resampling ◄─── 4. Saitou-Nei NJ
   (NWK + Dark Mode PNG)  (100 replicates, splits)      (Dynamic Q-Matrix)
```

### Supported Datasets (`--dataset`)
1. **`col1a1`**: 8 vertebrate comparative species across ~450 million years of evolutionary divergence (*Homo sapiens*, *Pan troglodytes*, *Mus musculus*, *Rattus norvegicus*, *Bos taurus*, *Gallus gallus*, *Danio rerio*, *Xenopus tropicalis*).
2. **`animals`**: 9 animal and outgroup taxa aligned across 839 positions from native MEGA format (`data/raw/animal_alignments/animals.meg`).
3. **`tro_seq`**: 30 zoonotic & vector epidemiological surveillance isolates (`data/raw/col1a1_sequences/TRO_Seq(30).fas`).

### Core Computational Algorithms

#### 1. Distance Calculation
- **JC69**: $d = -\frac{3}{4} \ln\left(1 - \frac{4}{3}p_{\text{total}}\right)$
- **K2P**: $d = -\frac{1}{2} \ln(1 - 2P - Q) - \frac{1}{4} \ln(1 - 2Q)$

#### 2. Neighbor-Joining ($Q$-Matrix Optimization)
- Iteratively evaluates:
  $$Q(i,j) = (n-2) D(i,j) - \sum_{k} D(i,k) - \sum_{k} D(j,k)$$
- Computes node branch lengths:
  $$d(i, u) = \frac{D(i,j)}{2} + \frac{\sum_k D(i,k) - \sum_k D(j,k)}{2(n-2)}, \quad d(j, u) = D(i,j) - d(i, u)$$

#### 3. Felsenstein (1985) Bootstrap Resampling (`bootstrap_neighbor_joining`)
- Resamples alignment columns with replacement ($B = 100$ replicates).
- Normalizes internal phylogenetic bipartitions into canonical unrooted splits:
  $$\text{canonical\_split}(C, S) = \begin{cases} S \setminus C & \text{if } s_0 \in C \\ C & \text{otherwise} \end{cases}$$
- Annotates bootstrap support frequencies directly into standard Newick trees and publication cladogram branches.

### Verification Suite (`tests/test_phylo.py`)
A comprehensive 11-test suite executed via `uv run pytest`:
- Analytical distance precision ($d=0.0$, expected K2P/JC69 values).
- Asymptotic saturation limits ($d \to \infty$ handling).
- Match, transition, transversion, and gap classification.
- Canonical split normalization invariance.
- Embedded bootstrap support within Newick formatting.
- Native `.meg` file parser extraction.
- Non-parametric bootstrap percentage bounded strictly in $[0, 100]\%$.

## Related Wiki Graph Nodes
- **Entities**: [[COL1A1]], [[MEGA12]], [[ITOL]], [[Homo_sapiens]], [[Pan_troglodytes]], [[Danio_rerio]], [[Phlebotomus_sergenti]]
- **Concepts**: [[Neighbor_Joining_Algorithm]], [[Pairwise_Evolutionary_Distance]], [[Bootstrap_Resampling_and_Phylogenetic_Support]], [[Phylogenetic_Reconstruction]], [[Substitution_Models]]
