import argparse as arg
from datetime import datetime
from database import *
import pybiomart as biomart
from sqlalchemy import create_engine
import pyranges
import pandas as pd
import itertools
from sqlalchemy import Table
from sqlalchemy import Column
from sqlalchemy import ForeignKey
from sqlalchemy import Integer
from sqlalchemy import String
from sqlalchemy import JSON
from sqlalchemy.orm import relationship
from sqlalchemy import MetaData
from sqlalchemy.orm import declarative_base



if __name__=="__main__":
    parser = arg.ArgumentParser(description='create a txdb-esque database, currently only creates a sqlite database')
    parser.add_argument('-g', '--gtf', help="gtf file", type=str, action="store")
    parser.add_argument('-o', '--output', help='path to the database')
    parser.add_argument('-m', '--biomart', help='biomart dataset name', default="hsapiens_gene_ensembl")
    parser.add_argument('-n', '--host', help="biomart host name", default="http://grch37.ensembl.org")
    parser.add_argument('-a', '--annotations', help="annotation columns from biomart comma separated text", default=None)
    args = parser.parse_args()

    gene_attr=["ensembl_gene_id", "external_gene_name",  "description", "gene_biotype"]
    tx_attr=["ensembl_transcript_id", "transcript_biotype"]

    dataset = biomart.Dataset(name=args.biomart,
                                  host=args.host)

    engine=create_engine("sqlite:///{}".format(args.output))
    genome_meta = MetaData(bind=engine)
    GenomeBase = declarative_base(metadata=genome_meta)


    class Chroms(GenomeBase):
        __tablename__ = "chrom"
        id = Column(String(100), primary_key=True)  # because of non cannonical and MT
        # I dont think I need anything else here


    class Genes(GenomeBase):
        __tablename__ = "genes"
        id = Column(String(100), primary_key=True, index=True)
        chrom = Column(String(), ForeignKey("chrom.id"))
        start = Column(Integer)
        end = Column(Integer)
        strand = Column(String(1))
        name = Column(String(100))
        description = Column(String())
        biotype = Column(String())
        user_annot = Column(JSON())
        chrom_rel = relationship("chrom", back_populates="genes", cascade="all")


    class Transcripts(GenomeBase):
        __tablename__ = "transcripts"
        id = Column(String(100), primary_key=True, index=True)
        gene = Column(String(), ForeignKey("genes.id"))
        start = Column(Integer)
        end = Column(Integer)
        biotype = Column(String)
        user_annot = Column(JSON())
        gene_rel = relationship("genes", back_populates="transcripts")
        exons_rel = relationship("exons", back_populates="transcripts")
        cds_rel = relationship("cdss", back_populates="transcripts")
        three_utr_rel = relationship("three_utr", back_populates="transcripts")
        five_utr_rel = relationship("five_utr", back_populates="transcripts")


    class Exons(GenomeBase):
        __tablename__ = "exons"
        id = Column(Integer, primary_key=True, index=True, autoincrement=True)
        exon_name = Column(String(100), index=True)
        transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
        start = Column(Integer)
        end = Column(Integer)
        rank = Column(Integer)
        transcript_rel = relationship("transcripts", back_populates="exons", cascade="all")


    class CDS(GenomeBase):
        __tablename__ = "cdss"
        id = Column(Integer, primary_key=True, index=True)
        transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
        exon_rank = Column(String(), ForeignKey("exons.id"))
        start = Column(Integer)
        end = Column(Integer)
        transcript_rel = relationship("transcripts", back_populates="cdss", cascade="all")
        exon_rel = relationship("exons", back_populates="cdss", cascade="all")


    class Three_UTR(GenomeBase):
        __tablename__ = "three_utr"
        id = Column(Integer, primary_key=True, index=True)
        transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
        start = Column(Integer)
        end = Column(Integer)
        transcript_rel = relationship("transcripts", back_populates="three_utr", cascade="all")
        exon_rel = relationship("exons", back_populates="three_utr", cascade="all")


    # should I create a single utr table and have a class?
    class Five_UTR(GenomeBase):
        __tablename__ = "five_utr"
        id = Column(Integer, primary_key=True)
        transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
        start = Column(Integer)
        end = Column(Integer)
        transcripts_rel = relationship("transcripts", back_populates="five_utr", cascade="all")
        exon_rel = relationship("exons", back_populates="five_utr", cascade="all")


    class Introns(GenomeBase):
        __tablename__ = "introns"
        id = Column(Integer, primary_key=True, index=True)
        transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
        start = Column(Integer)
        end = Column(Integer)
        transcripts_rel = relationship("transcripts", back_populates="introns", cascade="all")


    GenomeBase.metadata.create_all(engine)

    gene_features=dataset.query(attributes=gene_attr)
    tx_features=dataset.query(attributes=tx_attr)

    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "Reading gtf file " + args.gtf)
    gtf=pyranges.readers.read_gtf(args.gtf, full=False)

    # start inserting into database
    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "adding data to the database")
    chroms=pd.DataFrame({"id":gtf.Chromosome.unique()})
    chroms.to_sql("chrom", engine, if_exists="append", index=False)

    genes=gtf[(gtf.Feature=="gene")].df.loc[:, ["gene_id", "Chromosome", "Start", "End", "Strand"]]
    genes=genes.merge(gene_features, how="left", left_on="gene_id", right_on="Gene stable ID").\
        drop(columns=["Gene stable ID"])
    genes=genes.rename(columns={"gene_id":"id", "Chromosome":"chrom", "Start":"start", "End":"end",
                  "Strand":"strand", "Gene name":"name", "Gene description":"description",
                  "Gene type":"biotype"})
    genes.to_sql("genes", engine, if_exists="append", index=False)


    txs=gtf[(gtf.Feature=="transcript")].df.loc[:, ["Start", "End", "gene_id", "transcript_id"]]
    txs=txs.merge(tx_features, how="left", left_on="transcript_id", right_on="Transcript stable ID").\
        drop(columns=["Transcript stable ID"])
    txs=txs.rename(columns={"transcript_id":"id", "Start":"start", "End":"end","Transcript type":"biotype",
                            "gene_id":"gene"})
    txs.to_sql("transcripts", engine, if_exists="append", index=False)


    exons=gtf[(gtf.Feature=="exon")]
    if len(exons)>0:
        exons=exons.df.loc[:, ["Start", "End", "transcript_id", "exon_id", "exon_number"]]
        exons=exons.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript",
                                 "exon_number":"rank", "exon_id":"exon_name"})
        exons.to_sql("exons", engine, if_exists="append", index=False)


    cdss=gtf[(gtf.Feature=="CDS")]
    if len(cdss)>0:
        cdss=cdss.df.loc[:, ["Start", "End", "transcript_id", "exon_number"]]
        cdss=cdss.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript",
                                 "exon_number":"exon_rank"})
        cdss.to_sql("cdss", engine, if_exists="append", index=False)

    three_utr=gtf[(gtf.Feature=="three_prime_utr")]
    if len(three_utr)>0:
        three_utr=three_utr.df.loc[:, ["Start", "End", "transcript_id"]]
        three_utr=three_utr.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript"})
        three_utr.to_sql("three_utr", engine, if_exists="append", index=False)


    five_utr=gtf[(gtf.Feature=="five_prime_utr")]
    if len(five_utr)>0:
        five_utr=five_utr.df.loc[:, ["Start", "End", "transcript_id"]]
        five_utr=five_utr.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript"})
        five_utr.to_sql("five_utr", engine, if_exists="append", index=False)

    # this is different because there are no "introns" in the gtf we need to infer them from the locations of exons
    introns=gtf.features.introns(by="transcript")
    if len(introns)>0:
        introns=introns.df.loc[:, ["Start", "End", "transcript_id"]]
        introns=introns.rename(columns={"Start":"start", "End":"end", "transcript_id":"transcript"})
        introns.to_sql("introns", engine, if_exists="append", index=False)

    print("[" + datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] " + "adding_annotations to the database")
    if args.annotations is not None:
        annots=args.annotations.split(",")

        annots.insert(0, "ensembl_gene_id")
        types=list(itertools.repeat(String, len(annots)))
        idx=[True] + list(itertools.repeat(False, len(annots)-1))
        annots_table=Table("gene_annot", genome_meta,
                           *(Column(col, type, index=idx)
                             for col, type, ix in zip(annots, types, idx)))

        annots_table.create()
        gene_annotations = dataset.query(attributes=annots)
        cols=gene_annotations.columns
        new_cols={}
        for key, value in zip(cols, annots):
            new_cols[key]=value
        gene_annotations=gene_annotations.rename(columns=new_cols)
        gene_annotations.to_sql("gene_annot", engine, if_exists="append", index=False)

