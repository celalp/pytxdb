# python version of txdb database

This is my humble attempt to create something similar to R's GenomicFeatures library. I'm trying to create a sqlite database for low memory footprint given a gtf. 
The databse contains several tables that map genomic features to different tables. I'm also querying biomart for additional columns. There is also a gene
level annotations table that is dynamically created by again querying biomary for arbitrary columns. The user is responsible for not picking too many columns otherwise
biomart might throw an error. I will fix that limitation later on. 

## Creating the database

```
python3 reate_database.py  -g {path to gtf} -o {path to database} -m {ensembl name}   -n http://www.ensembl.org -a {comma,separated,ensembl_fields}
```

This creates a sqlite databse with the following schema:

Keep in mind that the columns in the gene_annot table will be dynamically created to match what's been provided in the `-a` field

## Querying the databse

After creating the database you would need to create a genome object

```python
from sqlalchemy import create_engine
import pybiomart as biomart
from pytxdb import *


db=create_engine("sqlite:///genome.db")
mart=biomart.Dataset(name={ensembl_name},
                    host={ensembl_url})
fasta="fasta_file_path" # chromosome names must match the gtf and be indexed
genome=Genome(db, fasta, mart)
```

There are several methods for getting different sequences these are

```python
Genome.genes()
Genome.transcript()
Genome.exons()
Genome.introns()
Genome.three_utr()
Genome.five_utr()
Genome.cds()
```

For genes (transcripts) you can specify gene (gene or transcript) names to return coordinates for
those only. The names needs to be a list even if you have only one thing. That will be addressed in the future

For other methods you can specify gene or transcript names but also need to tell whether you are 
looking for a gene or a transcript using the `level` argument. 

In the genes table I'm automatically querying gene symbol, description biotype. Additionally using the
other annotations that are requested in the gene_annot table you can search for ensembl ids with the
`search_gene`. You can search all fields or specify a column. This can be exact text or regex if `regex` parameter
is set to `True`. 


Each of the genomic feature queries (genes, exons etc.) either returns a pandas `DataFrame` or pyranges
`pyranges`. The latter can be used to return sequences from the fasta file using the `get_sequence` function. 
The returned sequences will be strand specific and biopython `Seq` object, if you provide more than 
one item in the genomic ranges you have the option to concat them (like coding sequences of a transcript to get the whole thing). 
You can also return amino acid sequences. I will implement a feature to translate starting
from the first `AUG` in a later version. Currently, you cannot specify an alternative codon table. That will be
addressed in the future as well. 

Additionally if you specify a pyranges object and a transcript name you can search for protein domains 
that match to that region. You will need a biomart connection and a valid query name like `pfam` using the
`domains` function. 

Finally, you can add arbitrary `dict` annotations to each gene or transcript these are stored as `JSON` in the
sqlite table but there isn't really a way to search them yet. That will be in the next release. 

Please let me know if there are any issues in the issues section 


