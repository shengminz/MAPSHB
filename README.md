# MAPSHB
The Machine Learning Assisted Prediction of Short Hydrogen Bonds (MAPSHB)

MAPSHB allows effective predictions of SHBs from a given protein structure. The MAPSHB model is designed for SHBs that are formed between two amino acids and have the donor residue on the side chain of an amino acid. This repository is a supplementary material of the paper "<a href="https://www.nature.com/articles/s41598-021-04306-4">Effective prediction of short hydrogen bonds in proteins via machine learning method</a>" and we encourage you to read it before using this pipeline.

**Colab notebook** <a href="https://colab.research.google.com/drive/1notF8VnttWgMMIkHiVNoWgEhQ9PjGkrK?authuser=1">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a> - `Using this Colab notebook to run MAPSHB prediction`

### Installation
1. Clone the repo
```bash
git clone git@github.com:JeffSHF/ColabDock.git
```
2. Install AmberTools23  
<b>Please refer to [JAX github](https://github.com/google/jax) page to install package corresponding to your CUDA and cudnn version.</b>  
Example:
```bash
# install jaxlib
pip install https://storage.googleapis.com/jax-releases/cuda11/jaxlib-0.3.8+cuda11.cudnn805-cp38-none-manylinux2014_x86_64.whl
# install jax
pip install jax==0.3.8
```
