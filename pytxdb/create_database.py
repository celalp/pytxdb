import argparse as arg
import os
from datetime import datetime

import pandas as pd
import pybiomart as biomart
import pyranges
import yaml

import utils
from database import *

if __name__ == "__main__":
    parser = arg.ArgumentParser(description='create a txdb-esque database, currently only creates a sqlite database')
    parser.add_argument('-g', '--gtf', help="gtf file", type=str, action="store")
    parser.add_argument('-y', '--yaml', help="config yaml file", type=str, action="store")
    args = parser.parse_args()

    if args.yaml is None:
        raise ValueError("no yaml file provided")

    if not os.path.isfile(args.yaml):
        raise FileNotFoundError("could not find {}".format(args.yaml))

    if args.gtf is None:
        raise ValueError("no no file provided")

    if not os.path.isfile(args.gtf):
        raise FileNotFoundError("could not find {}".format(args.gtf))

    with open(args.yaml) as y:
        params = yaml.safe_load(y)

    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Creating database " +
          params["output"]["name"])
    engine = utils.dict_to_engine(params["output"])
    genome_meta = MetaData(bind=engine)
    GenomeBase.metadata.create_all(engine)

    genome_meta.reflect()
    gene_annots_table = utils.dict_to_table(params["gene_annotations"], "gene_annotations", genome_meta)
    tx_annots_table = utils.dict_to_table(params["transcript_annotations"], "transcript_annotations", genome_meta)

    gene_annots_table.create()
    tx_annots_table.create()

    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Reading gtf file " + args.gtf)
    gtf = pyranges.readers.read_gtf(args.gtf, full=False)

    # start inserting into database
    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Adding data to the database")
    chroms = pd.DataFrame({"id": gtf.Chromosome.unique()})
    chroms.to_sql("chrom", engine, if_exists="append", index=False)

    genes = gtf[(gtf.Feature == "gene")].df.loc[:, ["gene_id", "Chromosome", "Start", "End", "Strand"]]
    genes = genes.rename(columns={"gene_id": "id", "Chromosome": "chrom", "Start": "start", "End": "end",
                                  "Strand": "strand"})
    genes.to_sql("genes", engine, if_exists="append", index=False)

    txs = gtf[(gtf.Feature == "transcript")].df.loc[:, ["Start", "End", "gene_id", "transcript_id"]]
    txs = txs.rename(columns={"transcript_id": "id", "Start": "start", "End": "end", "gene_id": "gene"})
    txs.to_sql("transcripts", engine, if_exists="append", index=False)

    exons = gtf[(gtf.Feature == "exon")]
    if len(exons) > 0:
        exons = exons.df.loc[:, ["Start", "End", "transcript_id", "exon_id", "exon_number"]]
        exons = exons.rename(columns={"Start": "start", "End": "end", "transcript_id": "transcript",
                                      "exon_number": "rank", "exon_id": "exon_name"})
        exons.to_sql("exons", engine, if_exists="append", index=False)

    cdss = gtf[(gtf.Feature == "CDS")]
    if len(cdss) > 0:
        cdss = cdss.df.loc[:, ["Start", "End", "transcript_id", "exon_number"]]
        cdss = cdss.rename(columns={"Start": "start", "End": "end", "transcript_id": "transcript",
                                    "exon_number": "exon_rank"})
        cdss.to_sql("cdss", engine, if_exists="append", index=False)

    three_utr = gtf[(gtf.Feature == "three_prime_utr")]
    if len(three_utr) > 0:
        three_utr = three_utr.df.loc[:, ["Start", "End", "transcript_id"]]
        three_utr = three_utr.rename(columns={"Start": "start", "End": "end", "transcript_id": "transcript"})
        three_utr.to_sql("three_utr", engine, if_exists="append", index=False)

    five_utr = gtf[(gtf.Feature == "five_prime_utr")]
    if len(five_utr) > 0:
        five_utr = five_utr.df.loc[:, ["Start", "End", "transcript_id"]]
        five_utr = five_utr.rename(columns={"Start": "start", "End": "end", "transcript_id": "transcript"})
        five_utr.to_sql("five_utr", engine, if_exists="append", index=False)

    # this is different because there are no "introns" in the gtf we need to infer them from the locations of exons
    introns = gtf.features.introns(by="transcript")
    if len(introns) > 0:
        introns = introns.df.loc[:, ["Start", "End", "transcript_id"]]
        introns = introns.rename(columns={"Start": "start", "End": "end", "transcript_id": "transcript"})
        introns.to_sql("introns", engine, if_exists="append", index=False)

    dataset = biomart.Dataset(name=params["dataset"], host=params["host"])

    if "gene_annotations" in params.keys():
        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Adding gene annotations from ensembl")
        gene_annotations = params["gene_annotations"]
        gene_annots = utils.mart_download(dataset, gene_annotations.keys(),
                                          common="ensembl_gene_id", error=params["error"])
        gene_annots.to_sql("gene_annotations", engine, index=False, if_exists="append")

    if "transcript_annotations" in params.keys():
        print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Adding transcript from ensembl")
        tx_annotations = params["transcript_annotations"]
        tx_annots = utils.mart_download(dataset, tx_annotations.keys(),
                                        common="ensembl_transcript_id", error=params["error"])

        tx_annots.to_sql("transcript_annotations", engine, index=False, if_exists="append")

    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Done")