---
title: Maximum Likelihood Phylogenetic Reconstruction and Model Selection
type: concept
tags: [maximum-likelihood, felsenstein, pruning-algorithm, bic, aic, modeltest]
created: 2026-10-27
updated: 2026-10-27
sources: ["[[COL1A1_Vertebrate_Phylogenetics_Dataset]]", "[[MEGA12]]", "[[Substitution_Models]]"]
---

# 📐 Maximum Likelihood Phylogenetic Reconstruction and Model Selection

## 📌 Theoretical Foundation

**Maximum Likelihood (ML)** in phylogenetics seeks the tree topology $\tau$ and branch lengths $\mathbf{t}$ that maximize the probability of observing the sequence alignment $\mathbf{D}$ given an explicit nucleotide substitution model $\mathbf{M}$:

$$\mathcal{L}(\tau, \mathbf{t}, \mathbf{M} \mid \mathbf{D}) = P(\mathbf{D} \mid \tau, \mathbf{t}, \mathbf{M})$$

Unlike distance-based methods like Neighbor-Joining ([[Neighbor_Joining_Algorithm]]), ML explicitly models every character transition at every site and naturally accommodates rate heterogeneity across lineages and sites.

---

## ⚙️ 1. Felsenstein's Pruning Algorithm (1981)

Directly calculating likelihood over all possible unobserved ancestral sequences requires enumerating $4^{M-2}$ ancestral state combinations per site. Felsenstein solved this exponential complexity through **post-order tree traversal (dynamic programming)**:

```
                  (Root Node 0)
                     /      \
               (Node 1)    (Node 2)
                 /   \       /   \
                A     B     C     D   (Observed Terminal Taxa)
```

For each site $u$ and internal node $k$ with child nodes $u_1$ and $u_2$:
$$L_k(s) = \left[ \sum_{s_1 \in \{A,C,G,T\}} P_{s \to s_1}(t_{k, u_1}) L_{u_1}(s_1) \right] \cdot \left[ \sum_{s_2 \in \{A,C,G,T\}} P_{s \to s_2}(t_{k, u_2}) L_{u_2}(s_2) \right]$$

Where:
* $L_k(s)$: Conditional likelihood that ancestral node $k$ had nucleotide state $s$.
* $P_{s \to s_1}(t)$: Transition probability over branch length $t$, defined by matrix exponential $\exp(\mathbf{Q} \cdot t)$ ([[Substitution_Models]]).

At the root node, the total likelihood for site $u$ is the weighted sum across base frequencies:
$$L^{(u)} = \sum_{s \in \{A,C,G,T\}} \pi_s L_{\text{root}}(s)$$

Total log-likelihood across all $N$ independent alignment sites:
$$\ln \mathcal{L} = \sum_{u=1}^{N} \ln L^{(u)}$$

---

## 🔬 2. Model Selection: BIC vs AIC Criterion

Before searching for the optimal ML tree, the best-fit substitution model must be selected to balance parameter complexity against likelihood gain:

$$\text{AIC} = 2k - 2 \ln \hat{\mathcal{L}}$$
$$\text{BIC} = k \ln(N) - 2 \ln \hat{\mathcal{L}}$$

Where:
* $k$: Number of free parameters in the model (e.g., base frequencies, transition/transversion rates, gamma shape $\alpha$, invariant proportion $I$).
* $N$: Total number of alignment sites.
* $\hat{\mathcal{L}}$: Maximum likelihood score achieved by the model.

> **Why BIC is preferred in MEGA 12:** Because $N$ is typically large in genomic datasets ($N > 5,000\text{ bp}$ in COL1A1), the $k \ln(N)$ penalty term in BIC is much stronger than $2k$ in AIC. This strictly prevents parameter overfitting while correctly identifying complex models (like **GTR+G+I**) only when genuinely supported by the empirical data.

---

## 🔗 Related Concepts & Sources
- [[Substitution_Models]] — Mathematical formulation of GTR, HKY, K2P, and JC69 rate matrices.
- [[Selection_Pressure_and_Rate_Heterogeneity]] — Empirical BIC parameter evaluation for COL1A1.
- [[Bootstrap_Resampling_and_Phylogenetic_Support]] — Testing ML tree branch support.
- [[MEGA12]] — Native implementation of ModelTest and heuristic tree search.
