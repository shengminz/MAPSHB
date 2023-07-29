# MAPSHB
The Machine Learning Assisted Prediction of Short Hydrogen Bonds (MAPSHB)

MAPSHB allows effective predictions of SHBs from a given protein structure. The MAPSHB model is designed for SHBs that are formed between two amino acids and have the donor residue on the side chain of an amino acid. This repository is a supplementary material of the paper "<a href="https://www.nature.com/articles/s41598-021-04306-4">Effective prediction of short hydrogen bonds in proteins via machine learning method</a>" and we encourage you to read it before using this pipeline.

**Colab notebook** <a href="https://colab.research.google.com/drive/1notF8VnttWgMMIkHiVNoWgEhQ9PjGkrK?authuser=1">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a> - `Using this Colab notebook to run MAPSHB prediction`

### Installation
1. Clone the repo
```bash
git clone https://github.com/shengminz/MAPSHB.git
```
2. Install AmberTools23  
```bash
conda create --name AmberTools23
conda activate AmberTools23
conda install -c conda-forge ambertools=23
```
3. Install R
```bash
conda install -c r r
conda install -c conda-forge r-gbm
```
4. Install Fortran
```bash
conda install -c conda-forge gfortran
```

