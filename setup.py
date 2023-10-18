import setuptools

setuptools.setup(
    name="example-app",
    version="0.0.1",
    install_requires=[
        "numpy",
        "rdkit",
        ],
    package_dir={"": "src"},
    packages=["biosynfoni"],
    python_requires=">=3.9",
    entry_points={"console_scripts": ["greeter = cli:main"]}
)
