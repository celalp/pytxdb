from setuptools import setup

setup(
    name='pytxdb',
    version='0.1.0',
    description="Python module for creating a sqlite database and searching for genomes inspired by R's GenomicFeatures and txdb",
    author='Alper Celik',
    author_email='alper.celik@sickkids.ca',
    packages=['pytxdb'],
    install_requires=["pandas","pyranges","pysam","biopython","sqlalchemy","numpy","pybiomart"],
    zip_safe=False,
    include_package_data=True
)