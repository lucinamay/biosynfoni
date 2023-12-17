<img width="800" alt="스크린샷 2023-10-19 오후 7 59 27" src="https://github.com/lucinamay/biosynfoni/assets/119406697/c2b32601-8a00-4520-b027-101206becf81">\
<span style="color:green"> 🌿 *a biosynformatic molecular fingerprint tailored to natural product chem- and bioinformatic research* 🌿</span>
\________________________________________________________________________________________


  **bi·o·syn·for·ma·tic**\
  /ˌbaɪ  oʊ  sɪn  fərˈ mæt ɪk/\
  *adjective Computers, Biochemistry*

  relating to biosynthetic information and biochemical logic.\
  as a concatenation of  *biosynthetic* and *bioinformatics*, it was coined\
  during the creation of `BioSynFoni`.

\_________________________________________________________________________________________


### usage 🌿

install the package by downloading the repository and in it run `pip install .` (requires `python > 3.9 `)\
(you can do this inside your desired conda environment if you want)\
then it's ready to use!\
you can call the biosynfoni by name as a stand-alone command-line program or import it as a package in your python scripts.\
most basic command line usage example:
- `biosynfoni "<SMILES>"` -> prints default biosynfoni version of the molecule
- `biosynfoni "<InChI>"` -> prints default biosynfoni version of the molecule
- `biosynfoni <molecule_supplier.sdf>` -> writes the biosynfonies of all the supplier's molecule to `csv` file

to explore other options, type `biosynfoni -h` or `biosynfoni --help`


\* for imports into jupyter notebook from a conda environment, you might need to additionally run\
  `%pip install biosynfoni` in the notebook before it can import biosynfoni.\
  if it does not work like that either, you might want to check your notebook's `sys.path` and\
   add `<path>/<to>/condaenv/<path>/<to>/<site-packages>`


### biosynfoni versions / options / substructure descriptions
current 'default version' == 'strictphenyl_1016' | no amino acids, no >6C rings\
default setting for blocking : block (change with flags: see `biosynfoni -h` or `biosynfoni --help` for more info\
newer versions not yet set as default version\
rest of info: under construction...\


### overview of package modules 🌿

converture | (c)O(n)VERTure | converts file of `InChI`'s or `SMILES` into a `RDKit Chem.Mol` object `sdf` file\
concerto_fp  |  Convert Chemical E-Representation TO fingerprint | converts `RDKit Chem.Mol` objects into their respective biosynfoni fingerprints\
def_biosynfoni 	| contains the definitions of biosynfoni (versions of each substructure key, \
                	versions of collections of substructure keys) as dictionaries\
inoutput 	|	collection of input and output handling functions that are reused between codes\


other files (e.g. in 'experiments') were used in collection and curation of data, and are written for specific scopes.\
reuse at own discretion.

### For feedback:
We welcome your feedback. For a streamlined feedback process, please open a new issue on github.com/biosynfoni if you cannot find your issue among existing ones. Please also make sure to add the appropriate subject in the title for clarity as follows, for easier understanding and faster help:

- **SMARTS representation** of a specific existing substructure
    - [SMARTS] `<issue title>`
- **adding/removing/merging substructures** within current biosynfoni 
    - [SUBSTRUCTURE] `<issue title>`
- **new substructure collection** for specific application 
    - [COLLECTION] `<issue title>`
- **detection algorithm** of the substructures within the molecule 
    - [DETECTION] `<issue title>`
- **output options** 
    - [OUTPUT] `<issue title>`
- **other feedback**  
    - [<1-keyword summary of feedback topic>] `<issue title>`

In addition, please explain any reasoning behind your feedback to help us speed up decisions on if and how to implement changes within the biosynfoni tool.

### FAIR principle-advised information
findability & accessability:\
-URL to the repository/code: https://github.com/lucinamay/biosynfoni\
interoperability:\
-dependencies: rdkit, numpy,tqdm, python>=3.9\ 
-related data: 
	- natural product database https://coconut.naturalproducts.net | https://doi.org/10.1186/s13321-020-00478-9, 
	-synthetic products were filtered out from zinc.: https://zinc.docking.org | https://pubs.acs.org/doi/10.1021/acs.jcim.0c00675 \
-related software: n.a.\
reusability\
-licence: MIT (see LICENCE file)\

badge of FAIRness:\
[![FAIR checklist badge](https://fairsoftwarechecklist.net/badge.svg)](https://fairsoftwarechecklist.net/v0.2?f=20&a=30112&i=20122&r=123)




