<img width="800" alt="스크린샷 2023-10-19 오후 7 59 27" src="https://github.com/lucinamay/biosynfoni/assets/119406697/c2b32601-8a00-4520-b027-101206becf81">\
<span style="color:green"> 🌿 *a biosynformatic molecular fingerprint tailored to natural product chem- and bioinformatic research* 🌿</span>


<p align="center">
    <a href="https://github.com/lucinamay/biosynfoni/actions/workflows/test-biosynfoni.yml">
        <img alt="Tests" src="https://github.com/lucinamay/biosynfoni/actions/workflows/test-biosynfoni.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/biosynfoni">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/biosynfoni" /></a>
    <a href="https://pypi.org/project/biosynfoni">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/biosynfoni" /></a>
    <a href="https://github.com/lucinamay/biosynfoni/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/cinemol" /></a>
    <a href='https://github.com/psf/black'>
        <img src='https://img.shields.io/badge/code%20style-black-000000.svg' alt='Code style: black' /></a>
    <!-- <a href="https://doi.org/10.5281/zenodo.11242217">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.11242217.svg" alt="DOI"></a> -->
    <a href="https://fairsoftwarechecklist.net/v0.2?f=20&a=30112&i=20122&r=123">
        <img src="https://fairsoftwarechecklist.net/badge.svg" alt="FAIR checklist badge"></a>
</p>

\________________________________________________________________________________________


  **bi·o·syn·for·ma·tic**\
  /ˌbaɪ  oʊ  sɪn  fərˈ mæt ɪk/\
  *adjective Computers, Biochemistry*

  relating to biosynthetic information and biochemical logic.\
  as a concatenation of  *biosynthetic* and *bioinformatics*, it was coined\
  during the creation of `BioSynFoni`.

\_________________________________________________________________________________________


### Getting started 🌿

#### Predict biosynthetic class

We have trained a biosynthetic class predictor on `biosynfoni` fingerprints. 

You can try out the predictor on your own molecules [here](https://moltools.bioinformatics.nl/biosynfoni)!

#### Installation

Biosynfoni requires Python 3.9 or later. RDKit is installed as a dependency when installing Biosynfoni.

To install the package, you can use pip:

```bash
pip install biosynfoni
```

Now you can import the `biosynfoni` package in your Python code or use the command line tool.

#### Usage in Python

Convert a SMILES string to a fingerprint:

```python
from biosynfoni import Biosynfoni
from rdkit import Chem

smi = <SMILES>
mol = Chem.MolFromSmiles(smi)
fp = Biosynfoni(mol).fingerprint  # returns biosynfoni's count fingerprint of the molecule
```

#### Usage in the command line

Create a fingerprint from a SMILES string:

```bash 
biosynfoni <SMILES>
```

Create a fingerprint from an InChI string:

```bash
biosynfoni <InChI>
```

Write the fingerprints of all molecules in an SDF file to a CSV file:

```bash
biosynfoni <molecule_supplier.sdf>
```

<!-- ### Preprint

#### Citation

If you use `biosynfoni` in your research, please cite our [preprint](https://chemrxiv.org/engage/chemrxiv/public-dashboard):

```bibtex
@article{nollen2025biosynfoni,
  title={Biosynfoni: A Biosynthesis-informed and Interpretable Lightweight Molecular Fingerprint},
  author={Nollen, Lucina-May, Meijer, David, Sorokina, Maria, and Van der Hooft, Justin J. J.},
  journal={chemRxiv},
  year={2025}
}
``` -->

#### Data availability

We created several biosynthetic class predictors for our manuscript, which can be downloaded from Zenodo [here](https://zenodo.org/records/14791239).

We have used data from the [COCONUT](https://coconut.naturalproducts.net) natural product database ([DOI](https://doi.org/10.1186/s13321-020-00478-9)) and [ZINC](https://zinc.docking.org) compound database ([DOI](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00675)). The parsed data used for the analysis in our manuscript can be downloaded from Zenodo [here](https://zenodo.org/records/14791205). 



