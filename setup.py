import setuptools

setuptools.setup(
    name="biosynfoni",
    version="0.0.1",
    #only one point for rdkit, numpy (e.g. make own smiles-> mol)
    install_requires=[
        "numpy",
        "rdkit",
        "tqdm",
        ],
    package_dir={"": "src"},
    packages=["biosynfoni"],
    python_requires=">=3.9",
    entry_points={"console_scripts": ["biosynfoni = biosynfoni.concerto_fp:main"]}
)
