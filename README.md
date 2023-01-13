# MicroCellClust

This repository contains implementations of solvers for the MicroCellClust optimisation problem presented in [1]. It is designed to search for rare cell-clusters that present a highly specific gene expression in scRNA-seq data. It has been implemented in the Scala programming language. This repository also contains an interface to run MicroCellClust in R.

* Version 1 contains the original solver [1], usable on small and middle-scale data

* Version 2 implements a new solver designed to solve the MicroCellClust optimization problem in large-scale single-cell data [2].



### References

[1] A. Gerniers, O. Bricard and P. Dupont (2021). MicroCellClust: mining rare and highly specific subpopulations from single-cell expression data. *Bioinformatics*, 37(19), 3220-3227. https://doi.org/10.1093/bioinformatics/btab239

[2] A. Gerniers, P. Dupont (2022). MicroCellClust 2: a hybrid approach for multivariate rare cell mining in large-scale single-cell data. In *2022 IEEE International Conference on Bioinformatics and Biomedicine (BIBM)*, 148-153.

