---
title: Multiple Sequence Alignment Algorithms in Phylogenomics
type: concept
tags: [msa, alignment, muscle, progressive-alignment, gap-penalties, bioinformatics]
created: 2026-10-27
updated: 2026-10-27
sources: ["[[COL1A1_Vertebrate_Phylogenetics_Dataset]]", "[[MEGA12]]", "[[Phylo_Pipeline_Script]]"]
---

# 🧬 Multiple Sequence Alignment Algorithms in Phylogenomics

## 📌 Overview

**Multiple Sequence Alignment (MSA)** is the essential foundational step in molecular evolutionary genetics. It arranges biological sequences (DNA, RNA, or amino acids) into a rectangular matrix such that homologous nucleotide positions (derived from a single ancestral site) occupy the same column.

In this pipeline, sequences retrieved from NCBI are aligned using **MUSCLE v3.8** (*Multiple Sequence Comparison by Log-Expectation*, Edgar 2004) embedded directly within [[MEGA12]].

---

## ⚙️ 1. MUSCLE v3.8 Algorithmic Phases

MUSCLE achieves both higher accuracy and substantially faster throughput than classical ClustalW by implementing a three-stage progressive refinement workflow:

```
Stage 1: Draft Progressive Alignment
├── Pairwise k-mer counting (word distances)
├── Fast UPGMA guide tree construction
└── Progressive profile-profile alignment along guide tree

Stage 2: Guide Tree Refinement
├── Re-calculate pairwise distances using Kimura distance over Stage 1 alignment
├── Construct refined guide tree via UPGMA
└── Re-align along new guide tree (retains only if SP score improves)

Stage 3: Iterative Refinement
├── Systematically delete edges in the guide tree to partition MSA into two sub-profiles
├── Re-align sub-profiles using profile-profile DP (Log-Expectation scoring)
└── Retain alignment changes that maximize Sum-of-Pairs (SP) score
```

---

## 🔬 2. Gap Scoring & The Log-Expectation (LE) Objective Function

The Log-Expectation score evaluates profile-profile matching between two aligned sequence sets $X$ and $Y$ at column $i$:

$$\text{LE}(X_i, Y_j) = (1 - f_X^G)(1 - f_Y^G) \sum_{A} \sum_{B} f_i^A f_j^B \ln \left( \frac{P_{AB}}{P_A P_B} \right) + \text{Gap Penalty}$$

Where:
* $f_i^A$: Observed frequency of amino acid or nucleotide $A$ at column $i$ in profile $X$.
* $f_X^G$: Frequency of gap characters in profile $X$ at that position.
* $P_{AB}$: Joint probability of substitution between characters $A$ and $B$.
* $P_A, P_B$: Background frequencies.

### Handling Coding Exons in COL1A1
In protein-coding genes like [[COL1A1]], indels (insertions and deletions) must preserve the open reading frame. In the vertebrate alignment:
* **Default Gap Open Penalty:** `-400` (discourages introducing gap openings into conserved structural domains).
* **Default Gap Extension Penalty:** `0` (permits longer contiguous gaps when whole codon deletions occur).
* **Codon Integrity:** Alignments are inspected to verify that nucleotide gaps occur in multiples of 3, preventing artificial reading-frame disruptions across exon junctions.

---

## 🔗 Related Concepts & Sources
- [[Phylogenetic_Reconstruction]] — Sequential role of alignment before substitution model evaluation.
- [[COL1A1]] — Vertebrate sequences aligned via MUSCLE.
- [[MEGA_File_Formats]] — Export formats (`.meg`, `.mas`, `.fasta`) storing aligned blocks.
- [[Pairwise_Evolutionary_Distance]] — Distance calculations derived directly from aligned columns.
