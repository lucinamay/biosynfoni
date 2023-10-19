
![33](https://github.com/lucinamay/biosynfoni/assets/119406697/93498121-b298-4cc7-ab3a-a7b7ca54b665)

<span style="color:green">a *biosynformatic* molecular fingerprint tailored to natural product chem- and bioinformatic research</span>
\________________________________________________________________________________


  **bi·o·syn·for·ma·tic**\
  /ˌbaɪ  oʊ  sɪn  fərˈ mæt ɪk/\
  *adjective Computers, Biochemistry*

  describes methods of retrieval and analysis of (bio)chemical information of\
  substances produced by living organisms, using biosynthetic information.\
  as a concatenation of  *biosynthetic* and *bioinformatics*, it was coined for the creation of BioSynFoni

\________________________________________________________________________________
 




_____ usage ____________________________________________________________________

install the package by downloading the repository and in it run "pip install ."
then it's ready to use!
you can call the biosynfoni by name as a stand-alone command-line program or import it as a package in your python scripts.
most basic command line usage example:
- biosynfoni "\<SMILES\>" -> prints default biosynfoni version of the molecule
- biosynfoni "\<InChI\>" -> prints default biosynfoni version of the molecule
- biosynfoni <molecule_supplier.sdf> -> writes the biosynfonies of all the supplier's molecule to csv file

\________________________________________________________________________________


_____ overview of package files ____________________________________________________________________

converture.py | (c)O(n)VERTure | converts file of InChI's or SMILES into a RDKit Chem.Mol object sdf file\
concerto_fp.py 	--  Converconverts RDKit Chem.Mol objects into their respective Biosynfoni fingerprints\
def_biosynfoni 	--  contains the definitions of Biosynfoni (versions of each substructure key, \
                	versions of collections of substructure keys) as dictionaries\
inoutput.py 	--	collection of input and output handling functions that are reused between codes\


Other files (e.g. in 'experiments') were used in collection and curation of data, and are written specifically for their data formats. Reuse at own discretion.


\____________________________________________________________________________________________________


