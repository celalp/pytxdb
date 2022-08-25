import argparse as arg
from datetime import datetime
from database import *

parser = arg.ArgumentParser(description='create a txdb-esque database, currently only creates a sqlite database')
parser.add_argument('-g', '--gtf', "gtf file",
                    type=str, action="store")
parser.add_argument('-o', '--output', help='path to the database')
parser.add_argument('-m', '--biomart', help='biomart dataset name', default="hsapiens_gene_ensembl")
parser.add_argument('-h', '--host', help="biomart host name", default="http://grch37.ensembl.org")
args = parser.parse_args()

gene_attr=["ensembl_gene_id", "external_gene_name",  "description", "gene_biotype"]
tx_attr=["ensembl_transcript_id", "transcript_biotype"]

dataset = biomart.Dataset(name=args.biomart,
                              host=args.host)

engine=create_engine("sqlite:///{output}".format(args.output))
GenomeBase.metadata.create_all(engine)

gene_features=dataset.query(attributes=gene_attr)
tx_features=dataset.query(attributes=tx_attr)

print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Reading gtf file " + args.gtf)
gtf=pyranges.readers.read_gtf(args.gtf, full=False)

# start inserting into database
print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "adding data to the database")
chroms=pd.DataFrame({"id":gtf.Chromosome.unique()})
chroms.to_sql("chrom", genome_engine, if_exists="append", index=False)

genes=gtf[(gtf.Feature=="gene")].df.loc[:, ["gene_id", "Chromosome", "Start", "End", "Strand"]]
genes=genes.merge(gene_features, how="left", left_on="gene_id", right_on="Gene stable ID").\
    drop(columns=["Gene stable ID"])
genes=genes.rename(columns={"gene_id":"id", "Chromosome":"chrom", "Start":"start", "End":"end",
              "Strand":"strand", "Gene name":"name", "Gene description":"description",
              "Gene type":"biotype"})
genes.to_sql("genes", genome_engine, if_exists="append", index=False)


txs=gtf[(gtf.Feature=="transcript")].df.loc[:, ["Chromosome", "Start", "End", "Strand", "gene_id", "transcript_id"]]
txs=txs.merge(tx_features, how="left", left_on="transcript_id", right_on="Transcript stable ID").\
    drop(columns=["Transcript stable ID"])
txs=txs.rename(columns={"transcript_id":"id", "Start":"start", "End":"end","Transcript type":"biotype",
                        "gene_id":"gene"})
txs.to_sql("transcripts", genome_engine, if_exists="append", index=False)


exons=gtf[(gtf.Feature=="exon")].df.loc[:, ["Chromosome", "Start", "End", "Strand", "transcript_id", "exon_id", "exon_number"]]
exons=exons.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript",
                             "exon_number":"rank", "exon_id":"exon_name"})
exons.to_sql("exons", genome_engine, if_exists="append", index=False)


cdss=gtf[(gtf.Feature=="CDS")].df.loc[:, ["Chromosome", "Start", "End", "Strand", "transcript_id", "exon_number"]]
cdss=cdss.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript",
                             "exon_number":"exon_rank"})
cdss.to_sql("cdss", genome_engine, if_exists="append", index=False)

three_utr=gtf[(gtf.Feature=="three_prime_utr")].df.loc[:, ["Chromosome", "Start", "End", "Strand", "transcript_id"]]
three_utr=three_utr.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript"})
three_utr.to_sql("three_utr", genome_engine, if_exists="append", index=False)

five_utr=gtf[(gtf.Feature=="five_prime_utr")].df.loc[:, ["Chromosome", "Start", "End", "transcript_id"]]
five_utr=five_utr.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript"})
five_utr.to_sql("five_utr", genome_engine, if_exists="append", index=False)

# this is different because there are no "introns" in the gtf we need to infer them from the locations of exons
introns=gtf.features.introns(by="transcript").df.loc[:, ["Chromosome", "Start", "End", "transcript_id"]]
introns=introns.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript"})
introns.to_sql("introns", genome_engine, if_exists="append", index=False)

print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "adding_annotations to the database")
#TODO make this a more dynamic process, maybe a json column
annot_cols=["ensembl_gene_id", "hgnc_symbol", "ucsc", "wikigene_name", "external_synonym"]
gene_annotations=dataset.query(attributes=annot_cols)
gene_annotations=gene_annotations[gene_annotations["Gene stable ID"].isin(genes.id.to_list())]
gene_annotations=gene_annotations.rename(columns={"Gene stable ID":"ens_id", "HGNC symbol":"hgnc",
                                                   "UCSC Stable ID":"ucsc", "WikiGene name":"wikigene",
                                                  "Gene Synonym":"synonyms"})
gene_annotations.to_sql("gene_annot", genome_engine, if_exists="append", index=False)
