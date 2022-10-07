from setuptools import setup, find_packages

setup(
    name='pytxdb',
    version='0.1.0',
    description="Python module for creating a sqlite database and searching for genomes inspired by R's GenomicFeatures and txdb",
    author='Alper Celik',
    author_email='alper.celik@sickkids.ca',
    packages=find_packages(),
    install_requires=["pandas","pyranges","pysam","biopython","sqlalchemy","numpy","pybiomart"],
    zip_safe=False,
    scripts=["pytxdb/create_database.py"],
    include_package_data=True
)
