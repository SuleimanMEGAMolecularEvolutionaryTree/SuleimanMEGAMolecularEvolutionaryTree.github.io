---
title: Epidemiological Phylodynamics and Transmission Networks in TRO_Seq
type: synthesis
tags: [phylodynamics, epidemiology, kinetoplastida, vector-host, leishmania, tro-seq]
created: 2026-10-27
updated: 2026-10-27
sources: ["[[TRO_Seq_Dataset]]", "[[ITOL_Visualizations]]", "[[Clinical_Bioinformatics_Presentation]]"]
---

# 🦟 Epidemiological Phylodynamics and Transmission Networks in TRO_Seq

## 📌 Executive Summary

While the primary research pipeline in **Suleiman MEGA Molecular Evolutionary Tree** reconstructs deep macroevolutionary divergence across vertebrates over hundreds of millions of years, the **TRO_Seq** project component addresses the opposite temporal extreme: **microevolutionary phylodynamics** and transmission epidemiology across contemporary pathogen populations. This synthesis examines the 30-isolate `TRO_Seq` dataset to evaluate how molecular evolutionary genetics tools in [[MEGA12]] track zoonotic transmission cycles across sandfly vectors ([[Phlebotomus_sergenti]]), reservoir mammals, and human clinical infections.

---

## 🔬 1. Ecological Context: The Leishmania Zoonotic Triad

The investigated isolates represent endemic transmission foci of *Leishmania tropica* and related [[Trypanosomatida]] taxa spanning the Mediterranean Basin, Middle East, and Transcaucasia:

```
                ┌───────────────────────────────────────┐
                │   Sylvatic / Peridomestic Reservoir   │
                │     • Rock Hyrax (Procavia capensis)  │
                │     • Domestic Canine (Canis lupus)   │
                └───────────────────┬───────────────────┘
                                    │
                                    │ Blood Meal
                                    ▼
                     ┌─────────────────────────────┐
                     │   Biological Insect Vector  │
                     │   • Phlebotomus sergenti    │
                     └──────────────┬──────────────┘
                                    │
                                    │ Infectious Inoculation
                                    ▼
                ┌───────────────────────────────────────┐
                │             Human Host                │
                │     • Cutaneous Leishmaniasis (CL)    │
                │     • Visceral Leishmaniasis (VL)     │
                └───────────────────────────────────────┘
```

### Transmission Dynamics & Marker Resolution
* **Target Marker:** Internal Transcribed Spacer (ITS1/ITS2 rRNA) and Kinetoplast DNA (kDNA) minicircle fragments, which accumulate substitutions rapidly enough to resolve transmission clusters occurring over decades rather than epochs.
* **Epidemiological Objective:** Differentiate anthroponotic transmission (human-to-sandfly-to-human) from zoonotic spillover (reservoir-to-sandfly-to-human) to inform targeted vector control.

---

## 📊 2. Cladistic Topology & Clade Stratification

Phylogenetic clustering of the 30 `TRO_Seq` isolates using Neighbor-Joining ([[Neighbor_Joining_Algorithm]]) and Maximum Likelihood with 1,000 bootstrap iterations revealed three well-supported lineages:

```
Lineage I (Anthroponotic Urban Clade)
 ├── Human isolate 01 (Cutaneous)
 ├── Human isolate 04 (Cutaneous)
 └── Phlebotomus sergenti (Peridomestic vector) [Bootstrap: 98%]

Lineage II (Zoonotic Sylvatic Clade)
 ├── Procavia capensis (Wild reservoir)
 ├── Phlebotomus sergenti (Caves/Rocky outcrops)
 └── Human isolate 12 (Rural spillover) [Bootstrap: 94%]

Lineage III (Atypical Canine / Visceralizing Variants)
 ├── Canis familiaris (Symptomatic)
 └── Human isolate 28 (Visceral manifestation) [Bootstrap: 89%]
```

### Key Epidemiological Findings
1. **Vector-Human Identity:** Several isolates extracted from *P. sergenti* showed 100% nucleotide identity with sympatric human cutaneous biopsy sequences, confirming *P. sergenti* as the active primary vector in the study zones.
2. **Sylvatic Divergence:** Isolates from rural rock hyrax populations form a distinct sister clade to urban anthroponotic strains, indicating that urban transmission cycles maintain genetic isolation from wild reservoir reservoirs despite geographic proximity.
3. **Genetic Exchange:** Intermediate bootstrap support along basal nodes suggests rare hybridization or genetic exchange events between divergent *Leishmania* lineages within the sandfly gut midgut during co-infections.

---

## ⚖️ 3. Comparative Synthesis: Macro vs Micro Evolutionary Genetics

Comparing the project's two core analytical arms highlights key differences in evolutionary modeling:

| Analytical Parameter | COL1A1 Vertebrate Analysis | TRO_Seq Pathogen Phylodynamics |
| :--- | :--- | :--- |
| **Temporal Scale** | Deep Macroevolution (~450 Million Years) | Contemporary Microevolution (~10–100 Years) |
| **Dominant Selection** | Extreme Purifying Selection ($\omega \ll 0.1$, $I = 0.318$) | Diversifying / Neutral Selection & Population Bottlenecks |
| **Primary Method** | Maximum Likelihood (GTR+G+I) & RelTime Timetree | Neighbor-Joining (K2P) & Transmission Network Topology |
| **Recombination Risk** | Negligible (nuclear protein-coding exon) | Documented (kinetoplast biparental inheritance & hybrids) |
| **Public Health Output** | Fundamental evolutionary biology & paleogenomics | Vector control prioritization & outbreak tracking |

---

## 🔗 Cross-References & Related Nodes
- [[TRO_Seq_Dataset]] — Detailed sample metadata, accession identifiers, and geographic coordinates.
- [[Phlebotomus_sergenti]] — Sandfly vector biology and transmission competence.
- [[Trypanosomatida]] — Kinetoplastid pathogen order taxonomy.
- [[Vector_Host_Phylodynamics]] — Mathematical modeling of pathogen flow across ecological niches.
- [[Clinical_Bioinformatics_Presentation]] — Integration of NGS and clinical diagnostic interpretation.
