# Establishing-Different-genes-in-smoking-and-non-smoking-population
# Differential Gene Expression Analysis using F-Test

## Overview

This project implements a statistical framework for identifying **differentially expressed genes** from microarray gene expression data using the **F-test**. The analysis compares nested linear models to determine whether gene expression significantly differs across experimental conditions.

The implementation is developed entirely in **Python** using **NumPy**, **Pandas**, and **SciPy**, providing a transparent implementation of the underlying statistical methodology.

---

## Objectives

* Load and preprocess microarray gene expression data.
* Construct design matrices representing the null and alternative hypotheses.
* Compute F-statistics for every gene.
* Estimate p-values using the F-distribution.
* Identify statistically significant genes.
* Visualize the distribution of F-statistics and p-values.

---

## Mathematical Formulation

Suppose the expression values for a gene are represented as

```
X = [x₁, x₂, ..., x₄₈]
```

where each element corresponds to one biological sample.

Two linear models are considered.

### Null Model

The null hypothesis assumes no interaction between the experimental factors.

```
H₀ : X = Nβ + ε
```

where

* **N** is the reduced design matrix,
* **β** represents the regression coefficients,
* **ε** denotes random experimental error.

---

### Alternative Model

The alternative hypothesis allows all treatment effects.

```
H₁ : X = Dβ + ε
```

where **D** is the full design matrix.

Since

```
N ⊂ D
```

the two models are nested.

---

### Projection Matrices

The projection matrices are computed as

```
P_N = N(NᵀN)⁻¹Nᵀ

P_D = D(DᵀD)⁻¹Dᵀ
```

The residual projection matrices become

```
R_N = I − P_N

R_D = I − P_D
```

where **I** is the identity matrix.

---

### F-Statistic

For every gene,

```
RSS₀ = XᵀR_NX

RSS₁ = XᵀR_DX
```

The F-statistic is

```
F = ((RSS₀ / RSS₁) − 1)
    ×
    ((n − rank(D)) /
     (rank(D) − rank(N)))
```

where

* **n = 48** samples,
* **rank(D)** is the rank of the full model,
* **rank(N)** is the rank of the reduced model.

A larger F-statistic indicates stronger evidence against the null hypothesis.

---

### P-value

The significance level is computed from the F-distribution

```
p = 1 − Fcdf(F ; d₁, d₂)
```

where

```
d₁ = rank(D) − rank(N)

d₂ = n − rank(D)
```

Genes satisfying

```
p < 0.05
```

are considered significantly differentially expressed.

---

## Methodology

1. Load the microarray expression dataset.
2. Square the expression values as part of preprocessing.
3. Construct the null (**N**) and full (**D**) design matrices.
4. Compute projection matrices.
5. Calculate the residual sum of squares for every gene.
6. Compute F-statistics.
7. Estimate p-values using the F-distribution.
8. Select statistically significant genes.
9. Visualize the statistical distributions.

---

## Features

* Differential gene expression analysis
* Matrix-based F-statistic computation
* Nested linear model comparison
* Automatic p-value estimation
* Significant gene identification
* Histogram of p-values
* Scatter plot of F-statistics

---

## Technologies Used

* Python
* NumPy
* Pandas
* SciPy
* Matplotlib
* Seaborn

---

## Project Structure

```text
.
├── data/
│   └── Raw_Data_GeneSpring.txt
│
├── assignment_3.py
├── README.md
└── figures/
    ├── histogram_of_pvalues.png
    └── F_statistics.png
```

---

## Installation

Clone the repository

```bash
git clone https://github.com/yourusername/differential-gene-expression.git

cd differential-gene-expression
```

Install dependencies

```bash
pip install numpy pandas scipy matplotlib seaborn statsmodels
```

---

## Running the Project

Update the dataset path inside `assignment_3.py` and execute

```bash
python assignment_3.py
```

The program will

* Load the gene expression dataset
* Compute F-statistics
* Calculate p-values
* Identify significant genes
* Generate statistical plots

---

## Output

The analysis produces

* F-statistic for every gene
* P-value for every gene
* List of significant genes
* Histogram of p-values
* Scatter plot of F-statistics

---

## Future Improvements

* Benjamini–Hochberg false discovery rate correction
* Volcano plots
* Heatmaps of significant genes
* Principal Component Analysis (PCA)
* Differential expression using linear models (LIMMA)
* RNA-seq count-based analysis with DESeq2 or edgeR

---

## License

This project is released under the MIT License.

---

## Author

Developed as part of a computational biology and bioinformatics assignment to implement differential gene expression analysis using statistical hypothesis testing and linear model comparison from first principles in Python.
