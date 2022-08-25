from sqlalchemy import Column
from sqlalchemy import ForeignKey
from sqlalchemy import Integer
from sqlalchemy import String
from sqlalchemy import JSON
from sqlalchemy.orm import relationship
from sqlalchemy import MetaData
from sqlalchemy.orm import declarative_base

genome_meta = MetaData()
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


# TODO this needs to be more flexible
class Gene_Annot(GenomeBase):
    __tablename__ = "gene_annot"
    id = Column(Integer, primary_key=True)
    ens_id = Column(String(), ForeignKey("genes.id"), index=True)
    hgnc = Column(String(), index=True)
    ucsc = Column(String(), index=True)
    wikigene = Column(String(), index=True)
    synonyms = Column(String(), index=True)
