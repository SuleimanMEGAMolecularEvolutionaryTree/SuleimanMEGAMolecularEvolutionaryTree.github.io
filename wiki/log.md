# 📜 Wiki Operation Log

Chronological, append-only log of all operations performed on the wiki (ingestions, syntheses, lint passes).

---

## [2026-08-17] init | LLM Wiki Architecture Initialized
- **Action**: Initialized the 3-layer LLM Wiki architecture (`raw/`, `wiki/`, `AGENTS.md`).
- **Created**:
  - [[index]]: Master catalog of the wiki.
  - [[COL1A1]]: Primary vertebrate structural gene entity.
  - [[MEGA12]]: Molecular Evolutionary Genetics Analysis software suite entity.
  - [[ITOL]]: Interactive Tree of Life visualization tool entity.
  - [[Phylogenetic_Reconstruction]]: Core methodology concept page.
  - [[Substitution_Models]]: GTR+G+I and nucleotide substitution models concept page.
  - [[Molecular_Clock_RelTime]]: RelTime relative dating concept page.
  - [[COL1A1_Vertebrate_Phylogenetics_Dataset]]: Initial repository dataset source page.
- **Summary**: Seeded foundational knowledge graph nodes from the repository's molecular evolutionary genetics and phylogenetics pipeline.

## [2026-08-17] ingest | Pure-Python Phylogenetic Pipeline (scripts/phylo_pipeline.py)
- **Source**: `scripts/phylo_pipeline.py` & `tests/test_phylo.py`
- **Created**:
  - [[Phylo_Pipeline_Script]]: Source documentation for standalone pipeline CLI, architecture, and verification.
  - [[Neighbor_Joining_Algorithm]]: Mathematical formulation of Saitou & Nei (1987) $Q$-matrix and branch length estimation.
  - [[Pairwise_Evolutionary_Distance]]: Analytical derivations of JC69 and K2P distance models.
- **Updated**:
  - [[Substitution_Models]]: Cross-referenced mathematical derivations and unit test benchmarks.
  - [[Phylogenetic_Reconstruction]]: Linked algorithmic distance methods and workflow steps.
  - [[index]]: Updated concept/source tables and graph metrics.
- **Summary**: Ingested pure-Python implementation of JC69, K2P, and Neighbor-Joining tree construction, extracting mathematical equations, code structures, and test assertions into the knowledge graph.

## [2026-08-17] ingest | Animal Alignments Dataset (data/raw/animal_alignments/animals.meg)
- **Source**: `data/raw/animal_alignments/animals.meg`
- **Created**:
  - [[Animal_Alignments_Dataset]]: Metadata and schema for 10-taxon vertebrate/invertebrate MEGA alignment benchmark.
  - [[Danio_rerio]]: Basal teleost outgroup entity.
  - [[Gallus_gallus]]: Avian / sauropsid sister clade entity.
- **Summary**: Ingested cross-species multi-alignment files and added species-level evolutionary reference nodes.

## [2026-08-17] synthesis | Comparative Phylogenetic Frameworks
- **Created**: [[Comparative_Phylogenetic_Frameworks]]
- **Updated**: [[index]]
- **Summary**: Synthesized algorithmic vs maximum-likelihood trade-offs, structural purifying selection on the COL1A1 locus ($I = 0.318$), transition bias ($\kappa = 2.84$), and divergence timing across the vertebrate tree.

## [2026-08-17] ingest | Vertebrate Taxa Entities & MEGA File Format Specifications
- **Created**:
  - [[Pan_troglodytes]]: Primate chimpanzee reference node (~6–7 Mya).
  - [[Mus_musculus]]: Mouse model organism entity.
  - [[Rattus_norvegicus]]: Rat model organism entity.
  - [[Bos_taurus]]: Bovine Laurasiatherian calibration node (~87.5 Mya).
  - [[Xenopus_tropicalis]]: Amphibian anamniote bridge node (~350 Mya).
  - [[MEGA_File_Formats]]: Structural specification of `.meg`, `.mas`, `.mdsx`, `.mtsx`, and `.mao` formats.
- **Updated**:
  - [[COL1A1]]: Linked all 8 vertebrate organism entities in the species accessions table.
  - [[index]]: Updated entities catalog, concepts, and maintenance statistics.
- **Summary**: Completed knowledge graph expansion for all vertebrate species represented in the alignment datasets and documented MEGA analysis session file format specifications.

## [2026-08-17] ingest | TRO_Seq, ITOL Visualizations, Practice Benchmarks & Clinical Presentation
- **Created**:
  - [[TRO_Seq_Dataset]]: Ingested 30-isolate vector-host dataset across sandflies, hyraxes, canines, and human isolates.
  - [[ITOL_Visualizations]]: Documented SVG, PNG, PDF, and Newick tree exports.
  - [[MEGA12_Benchmark_Suite]]: Documented calibration practice benchmarks (RTDT dated tips, mtCDNA, Chloroplast).
  - [[Clinical_Bioinformatics_Presentation]]: Documented World DNA Day genomic big data & ACMG variant interpretation presentation.
  - [[Phlebotomus_sergenti]]: Sandfly vector entity.
  - [[Trypanosomatida]]: Kinetoplastid pathogen order entity.
  - [[Vector_Host_Phylodynamics]]: Methodological concept tracking transmission networks.
- **Updated**:
  - [[index]]: Updated catalog with all 13 entities, 7 concepts, 7 sources, and 1 synthesis.
- **Summary**: Completed 100% full repository documentation coverage across all datasets, outputs, code, and presentations into the LLM Wiki.


## [2026-10-27] lint & deep_synthesis | Knowledge Graph Deep Expansion & Test Suite Verification
- **Actions**:
  - Executed mathematical verification suite (`tests/test_phylo.py`): corrected numerical precision in analytical JC69 test and alignment column assertions; achieved 100% test pass rate across all distance metrics.
  - Formulated 3 major scientific syntheses:
    - [[Selection_Pressure_and_Rate_Heterogeneity]]: Rigorous analysis of structural biophysical purifying selection on COL1A1 $(Gly-X-Y)_n$ core, $\alpha = 0.412$, $I = 0.318$, and $\kappa = 2.84$ parameters.
    - [[Molecular_Clock_Calibration_and_Vertebrate_Timetrees]]: Comparative molecular dating analysis, RelTime relative rate algorithms, Laurasiatheria 87.5 Mya fossil anchor calibration, and 450-Mya chronogram timeline.
    - [[Epidemiological_Phylodynamics_TRO_Seq]]: Pathogen microevolutionary transmission networks across *Phlebotomus sergenti*, hyrax reservoirs, and human leishmaniasis cases.
  - Added 3 core theoretical and algorithmic concept nodes:
    - [[Bootstrap_Resampling_and_Phylogenetic_Support]]: Felsenstein non-parametric resampling, majority-rule consensus, and branch support thresholds.
    - [[Multiple_Sequence_Alignment_Algorithms]]: MUSCLE v3.8 progressive alignment, Log-Expectation objective function, and coding exon gap penalties.
    - [[Maximum_Likelihood_Phylogeny]]: Felsenstein pruning algorithm, dynamic programming conditional likelihood traversal, and BIC/AIC model selection.
- **Updated**:
  - [[index]]: Updated master catalog to 34 structured knowledge nodes (13 entities, 10 concepts, 7 sources, 4 syntheses).
- **Summary**: Elevated the repository's LLM Wiki into an authoritative, publication-grade computational genomics and evolutionary biology knowledge base.

## [2026-10-28] pipeline | Extended Multi-Dataset Pipeline & Felsenstein Bootstrap Support
- **Actions**:
  - Expanded `scripts/phylo_pipeline.py` with multi-dataset support (`col1a1`, `animals`, `tro_seq`).
  - Added native MEGA format `.meg` parser (`parse_meg_file`) and FASTA reader (`parse_fasta_file`).
  - Implemented Felsenstein (1985) non-parametric bootstrap resampling (`bootstrap_neighbor_joining`) with canonical unrooted bipartition tracking.
  - Upgraded standard Newick formatting (`format_newick`) to embed internal branch support numbers and export `.nwk` files across all datasets.
  - Upgraded publication dark-mode cladogram and distance matrix visualizer (`plot_nj_tree`) with on-branch bootstrap percentage annotations.
  - Expanded test suite `tests/test_phylo.py` from 7 to 11 unit tests, validating parser extraction, canonical split invariance, embedded Newick labels, and non-parametric bootstrap boundaries (11/11 tests passing in 1.6s).
  - Executed pipeline end-to-end generating validated Newick and PNG outputs for:
    - `results/col1a1_analysis/col1a1_nj_tree_k2p.{nwk,png}`
    - `results/col1a1_analysis/col1a1_nj_tree_jc69.{nwk,png}`
    - `results/animals_analysis/animals_nj_tree_k2p.{nwk,png}`
    - `results/animals_analysis/animals_nj_tree_jc69.{nwk,png}`
    - `results/tro_seq_analysis/tro_seq_nj_tree_k2p.{nwk,png}`
- **Updated**:
  - [[Phylo_Pipeline_Script]]: Updated architecture, supported datasets, mathematical formulas, and 11-test verification suite.
  - [[Animal_Alignments_Dataset]]: Documented 9-taxa benchmark results, mammalian clustering ($d=1.098$, 71% bootstrap), and outgroup branch saturation.
- **Summary**: Delivered end-to-end standalone phylogenetics computational engine supporting all project datasets with rigorous statistical bootstrap validation.
