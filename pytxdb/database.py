from sqlalchemy import Column, ForeignKey, Integer, String, JSON
from sqlalchemy.orm import declarative_base


GenomeBase = declarative_base()

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
    user_annot = Column(JSON())


class Transcripts(GenomeBase):
    __tablename__ = "transcripts"
    id = Column(String(100), primary_key=True, index=True)
    gene = Column(String(), ForeignKey("genes.id"))
    start = Column(Integer)
    end = Column(Integer)
    user_annot = Column(JSON())


class Exons(GenomeBase):
    __tablename__ = "exons"
    id = Column(Integer, primary_key=True, index=True, autoincrement=True)
    exon_name = Column(String(100), index=True)
    transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
    start = Column(Integer)
    end = Column(Integer)
    rank = Column(Integer)


class CDS(GenomeBase):
    __tablename__ = "cdss"
    id = Column(Integer, primary_key=True, index=True)
    transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
    exon_rank = Column(Integer)
    start = Column(Integer)
    end = Column(Integer)


class Three_UTR(GenomeBase):
    __tablename__ = "three_utr"
    id = Column(Integer, primary_key=True, index=True)
    transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
    start = Column(Integer)
    end = Column(Integer)


# should I create a single utr table and have a class?
class Five_UTR(GenomeBase):
    __tablename__ = "five_utr"
    id = Column(Integer, primary_key=True)
    transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
    start = Column(Integer)
    end = Column(Integer)


class Introns(GenomeBase):
    __tablename__ = "introns"
    id = Column(Integer, primary_key=True, index=True)
    transcript = Column(String(), ForeignKey("transcripts.id"), index=True)
    start = Column(Integer)
    end = Column(Integer)

