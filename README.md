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


\
_____ usage ____________________________________________________________________

install the package by downloading the repository and in it run `pip install .`\
then it's ready to use!\
you can call the biosynfoni by name as a stand-alone command-line program or import it as a package in your python scripts.\
most basic command line usage example:
- `biosynfoni "<SMILES>"` -> prints default biosynfoni version of the molecule
- `biosynfoni "<InChI>"` -> prints default biosynfoni version of the molecule
- `biosynfoni <molecule_supplier.sdf>` -> writes the biosynfonies of all the supplier's molecule to `csv` file

to explore other options, type `biosynfoni -h` or `biosynfoni --help`

\________________________________________________________________________________


_____ overview of package files ____________________________________________________________________

converture.py | (c)O(n)VERTure | converts file of `InChI`'s or `SMILES` into a `RDKit Chem.Mol` object `sdf` file\
concerto_fp.py  |  Convert Chemical E-Representation TO fingerprint | converts `RDKit Chem.Mol` objects into their respective biosynfoni fingerprints\
def_biosynfoni 	| contains the definitions of biosynfoni (versions of each substructure key, \
                	versions of collections of substructure keys) as dictionaries\
inoutput.py 	|	collection of input and output handling functions that are reused between codes\


other files (e.g. in 'experiments') were used in collection and curation of data, and are written specifically for their data formats. Reuse at own discretion.


\____________________________________________________________________________________________________


