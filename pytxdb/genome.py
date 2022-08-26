import os
import numpy as np
import pandas as pd
import pyranges
import pysam
import sqlalchemy as sql
from Bio.Seq import Seq
from sqlalchemy.orm import Session
from .utils import check_results

class Genome:
    def __init__(self, db, fasta=None, mart=None):
        """
        get the database connection, and fasta connection for getting sequences
        :param db: genome database, this is a connection
        :param fasta: fasta file for getting sequences optional.
        :param mart pybiomart connection optional
        """
        self.db = db
        if not os.path.isfile(fasta) and fasta is not None:
            raise FileNotFoundError("{} does not exist".format(fasta))
        else:
            self.fasta = pysam.Fastafile(fasta)
        self.mart = mart
        self.session = Session(self.db)
        self.metadata = sql.MetaData(self.db)
        self.metadata.reflect(bind=self.db)

    def search_gene(self, term, where=None, regex=False, return_fields=False):
        """
        searches for ensembl id based on some information like gene name can use regexes
        :param term: what to search for
        :param where: optional search field
        :param regex: is this a regex default False
        :return: returns a list of strings (ens ids)
        """

        if return_fields:
            print("Available fields are: name, description, biotype, hgnc, ucsc, wikigene, synonyms")
        else:

            gene_table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)
            annot_table = sql.Table("gene_annot", self.metadata, autoload=True, autoload_with=self.db)

            query1 = sql.select(gene_table.c.id, gene_table.c.name, gene_table.c.description,
                                gene_table.c.biotype)

            res1 = self.session.execute(query1).fetchall()
            df1 = pd.DataFrame(res1)
            columns = list(self.session.execute(query1).keys())
            df1.columns = columns

            query2 = sql.select(annot_table.c.ens_id, annot_table.c.hgnc, annot_table.c.ucsc,
                                annot_table.c.wikigene, annot_table.c.synonyms)
            res2 = self.session.execute(query2).fetchall()
            df2 = pd.DataFrame(res2)
            columns = list(self.session.execute(query2).keys())
            df2.columns = columns

            merged = df1.merge(df2, how="left", right_on="ens_id", left_on="id"). \
                drop(columns="ens_id")

            mask_array = np.zeros_like(merged, dtype=bool)
            if where is None:
                for i in range(len(merged.columns)):
                    mask_array[:, i] = merged.iloc[:, i]. \
                        str.contains(term, case=False, regex=regex, na=False)
                mask = np.sum(mask_array, axis=1) > 0

            else:
                mask = merged[where].str.contains(term, case=False, regex=regex, na=False)

            genes = merged["id"][mask].drop_duplicates()
            return genes

    def genes(self, names=None, df=False, additional_features=True):
        """
        get gene coordinates
        :param names: list of ens ids if none everything
        :param df: return a dataframe or pyranges
        :param additional_features: add additional features to the search results
        :return: returns a dataframe or pyranges object
        """
        table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)
        if names is None:
            query = sql.select(table)
        else:
            query = sql.select(table).where(table.c.id.in_(names))
        results = check_results(self.session.execute(query).fetchall())

        columns = list(self.session.execute(query).keys())

        dat = pd.DataFrame(results)
        dat.columns = columns

        if df:
            return dat
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=dat.iloc[:, 1],
                                            starts=dat.iloc[:, 2],
                                            ends=dat.iloc[:, 3],
                                            strands=dat.iloc[:, 4])
            if additional_features:
                gr.id = dat["id"]
                gr.name = dat["name"]
                gr.description = dat["description"]
                gr.biotype = dat["biotype"]
            return gr

    def transcripts(self, names=None, df=False, additional_features=True):
        """
        return transcripts
        :param names: list of ens ids if none everything
        :param df: return a df otherwise pyranges
        :param additional_features: add additional features to the search results
        :return: a df or pyranges
        """
        tx_table = sql.Table("transcripts", self.metadata, autoload=True, autoload_with=self.db)
        gene_table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)
        # need to get chrom and strand
        subq = sql.select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()

        query = sql.select(tx_table, subq.c.chrom, subq.c.strand). \
            join(subq, tx_table.c.gene == subq.c.id)

        if names is not None:
            query = query.where(tx_table.c.gene.in_(names))

        results = check_results(self.session.execute(query).fetchall())
        dat = pd.DataFrame(results)
        dat.columns = list(self.session.execute(query).keys())

        if df:
            return dat
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=dat.iloc[:, 6],
                                            starts=dat.iloc[:, 2],
                                            ends=dat.iloc[:, 3],
                                            strands=dat.iloc[:, 7])
            if additional_features:
                gr.tx_id = dat.id
                gr.gene_id = dat.gene
                gr.biotype = dat.biotype
            return gr

    def exons(self, names=None, level="transcript", df=False, additional_features=True):
        """
        return exons
        :param names names of transcripts (or genes)
        :param level for names either "gene" or "transcript"
        :param df: return a df?
        :param additional_features: add additional features to the search results
        :return: a df or pyranges
        """
        exon_table = sql.Table("exons", self.metadata, autoload=True, autoload_with=self.db)
        tx_table = sql.Table("transcripts", self.metadata, autoload=True, autoload_with=self.db)
        gene_table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)

        subq1 = sql.select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()
        subq2 = sql.select(tx_table, tx_table.c.gene, subq1.c.chrom, subq1.c.strand). \
            join(subq1, tx_table.c.gene == subq1.c.id).subquery()
        query = sql.select(exon_table, subq2.c.gene, subq2.c.chrom, subq2.c.strand). \
            join(subq2, exon_table.c.transcript == subq2.c.id)

        if names is not None:
            if level == "gene":
                query = query.where(subq2.c.gene.in_(names))
            elif level == "transcript":
                query = query.where(subq2.c.id.in_(names))
            else:
                raise ValueError("level can only be gene or transcript")

        results = check_results(self.session.execute(query).fetchall())
        dat = pd.DataFrame(results)
        dat.columns = list(self.session.execute(query).keys())
        dat = dat.drop(columns="id")

        if df:
            return dat
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=dat.iloc[:, 6],
                                            starts=dat.iloc[:, 2],
                                            ends=dat.iloc[:, 3],
                                            strands=dat.iloc[:, 7])
            if additional_features:
                gr.tx_id = dat.transcript
                gr.gene_id = dat.gene
                gr.rank = dat["rank"]
                gr.exon_name = dat.exon_name

            return gr

    def cds(self, names=None, df=False, level="gene", additional_features=True):
        """
        same as exons
        :param names:
        :param df:
        :param level
        :param additional_features
        :return:
        """
        cds_table = sql.Table("cdss", self.metadata, autoload=True, autoload_with=self.db)
        tx_table = sql.Table("transcripts", self.metadata, autoload=True, autoload_with=self.db)
        gene_table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)

        subq1 = sql.select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()
        subq2 = sql.select(tx_table, tx_table.c.gene, subq1.c.chrom, subq1.c.strand). \
            join(subq1, tx_table.c.gene == subq1.c.id).subquery()

        query = sql.select(cds_table, subq2.c.gene, subq2.c.chrom, subq2.c.strand). \
            join(subq2, cds_table.c.transcript == subq2.c.id)

        if names is not None:
            if level == "gene":
                query = query.where(subq2.c.gene.in_(names))
            elif level == "transcript":
                query = query.where(subq2.c.id.in_(names))
            else:
                raise ValueError("level can only be gene or transcript")

        results = check_results(self.session.execute(query).fetchall())
        dat = pd.DataFrame(results)
        dat.columns = list(self.session.execute(query).keys())
        dat = dat.drop(columns=["id", "exon_name"])

        if df:
            return dat
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=dat.iloc[:, 5],
                                            starts=dat.iloc[:, 2],
                                            ends=dat.iloc[:, 3],
                                            strands=dat.iloc[:, 6])
            if additional_features:
                gr.tx_id = dat.transcript
                gr.gene_id = dat.gene
                gr.exon_rank = dat.exon_rank

            return gr

    def introns(self, names=None, df=False, level="gene", additional_features=True):
        """
        same as exons
        :param names:
        :param df:
        :param level
        :return:
        """
        intron_table = sql.Table("introns", self.metadata, autoload=True, autoload_with=self.db)
        tx_table = sql.Table("transcripts", self.metadata, autoload=True, autoload_with=self.db)
        gene_table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)

        subq1 = sql.select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()
        subq2 = sql.select(tx_table, tx_table.c.gene, subq1.c.chrom, subq1.c.strand). \
            join(subq1, tx_table.c.gene == subq1.c.id).subquery()

        query = sql.select(intron_table, subq2.c.gene, subq2.c.chrom, subq2.c.strand). \
            join(subq2, intron_table.c.transcript == subq2.c.id)

        if names is not None:
            if level == "gene":
                query = query.where(subq2.c.gene.in_(names))
            elif level == "transcript":
                query = query.where(subq2.c.id.in_(names))
            else:
                raise ValueError("level can only be gene or transcript")

        results = self.session.execute(query).fetchall()
        dat = pd.DataFrame(results)
        dat.columns = list(self.session.execute(query).keys())
        dat = dat.drop(columns="id")

        if df:
            return dat
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=dat.iloc[:, 4],
                                            starts=dat.iloc[:, 1],
                                            ends=dat.iloc[:, 2],
                                            strands=dat.iloc[:, 5])
            if additional_features:
                gr.tx_id = dat.transcript
                gr.gene_id = dat.gene
            return gr

    def three_utr(self, names=None, df=False, level="gene", additional_features=True):
        """
        same as exons
        :param names:
        :param df:
        :return:
        """
        utr_table = sql.Table("three_utr", self.metadata, autoload=True, autoload_with=self.db)
        tx_table = sql.Table("transcripts", self.metadata, autoload=True, autoload_with=self.db)
        gene_table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)

        subq1 = sql.select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()
        subq2 = sql.select(tx_table, tx_table.c.gene, subq1.c.chrom, subq1.c.strand). \
            join(subq1, tx_table.c.gene == subq1.c.id).subquery()

        query = sql.select(utr_table, subq2.c.gene, subq2.c.chrom, subq2.c.strand). \
            join(subq2, utr_table.c.transcript == subq2.c.id)

        if names is not None:
            if level == "gene":
                query = query.where(subq2.c.gene.in_(names))
            elif level == "transcript":
                query = query.where(subq2.c.id.in_(names))
            else:
                raise ValueError("level can only be gene or transcript")

        results = self.session.execute(query).fetchall()
        dat = pd.DataFrame(results)
        dat.columns = list(self.session.execute(query).keys())

        if df:
            return dat
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=dat.iloc[:, 5],
                                            starts=dat.iloc[:, 2],
                                            ends=dat.iloc[:, 3],
                                            strands=dat.iloc[:, 6])
            if additional_features:
                gr.tx_id = dat.transcript
                gr.gene_id = dat.gene
            return gr

    def five_utr(self, names=None, df=False, level="gene", additional_features=True):
        """
        same as exons
        :param names:
        :param level:
        :param df:
        :return:
        """
        utr_table = sql.Table("five_utr", self.metadata, autoload=True, autoload_with=self.db)
        tx_table = sql.Table("transcripts", self.metadata, autoload=True, autoload_with=self.db)
        gene_table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)

        subq1 = sql.select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()
        subq2 = sql.select(tx_table, tx_table.c.gene, subq1.c.chrom, subq1.c.strand). \
            join(subq1, tx_table.c.gene == subq1.c.id).subquery()

        query = sql.select(utr_table, subq2.c.gene, subq2.c.chrom, subq2.c.strand). \
            join(subq2, utr_table.c.transcript == subq2.c.id)

        if names is not None:
            if level == "gene":
                query = query.where(subq2.c.gene.in_(names))
            elif level == "transcript":
                query = query.where(subq2.c.id.in_(names))
            else:
                raise ValueError("level can only be gene or transcript")

        results = self.session.execute(query).fetchall()
        dat = pd.DataFrame(results)
        dat.columns = list(self.session.execute(query).keys())

        if df:
            return dat
        else:
            gr = pyranges.pyranges.PyRanges(chromosomes=dat.iloc[:, 5],
                                            starts=dat.iloc[:, 2],
                                            ends=dat.iloc[:, 3],
                                            strands=dat.iloc[:, 6])
            if additional_features:
                gr.tx_id = dat.transcript
                gr.gene_id = dat.gene
            return gr

    def add_annotation(self, to, annot, gene=True):
        """
        add some annotations to the genes or transcripts table
        :param to: id of the gene or transcript
        :param annot: dict of annotations
        :param gene: is this a gene otherwise transcript
        :return: nothing adds the dict to the json field of either gene or transcript table
        """
        if gene:
            table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)
        else:
            table = sql.Table("transcripts", self.metadata, autoload=True, autoload_with=self.db)

        command = sql.insert(table).values(user_annot=annot).where(table.c.id == to)
        self.session.execute(command)
        self.session.commit()

    def get_sequence(self, gr, type="nuc", concat=True):
        """
        get the sequences from genomic ranges
        :param gr: a genomic ranges object
        :param type: "nuc" or "aa"
        :return: Seq object
        """
        if self.fasta is None:
            raise ValueError("Need to specify fasta location")
        seqs = []
        for i in range(len(gr)):
            sequence = self.fasta.fetch(reference=gr.Chromosome.to_list()[i],
                                        start=gr.Start.to_list()[i],
                                        end=gr.End.to_list()[i])
            sequence = Seq(sequence)
            if gr.Strand[i] == "+":
                seqs.append(sequence)
            else:
                sequence = sequence.reverse_complement()
                seqs.append(sequence)

        if type == "aa":
            seqs = [seq + ("N" * (3 - (len(seq) % 3))) for seq in seqs if len(seq) % 3 != 0]
            seqs = [seq.translate() for seq in seqs]

        if concat and len(seqs) > 1:
            seqs = seqs[0].join([seqs[i] for i in range(1, len(seqs))])

        return seqs

    def domains(self, gr, feature_name, tx):
        """
       if a given genomic ranges falls into a gene (or multiple genes) finds the regions that are covered and extracts
       proteins domains that are covered in that region
        :param gr: a pyranges to search the genome
        :param feature_name, what kind of feature like pfam domain or some other thing, uses biomart
        :param tx: transcript name, query only one transcript and feature at a time so ensembly servers
        respond in a reasonable amount of time
        :return:
        """
        features = ["gene3d", "hamap", "hmmpanther", "pfam", "pirsf", "prints", "scanprosite", "pfscan",
                    "smart", "superfamily", "tigrfam", "interpro", "ncoils", "seg", "signalp", "tmhmm"]
        if self.mart is None:
            raise ValueError("Need to specify mart connection")

        if feature_name not in features:
            raise ValueError("features can be one of the following: {}".format(",".join(features)))

        coding_regions = self.cds(names=[tx], level="transcript", df=False)
        strand = coding_regions.Strand.unique()[0]

        if strand != gr.Strand:
            print("region provided and the transcript are on different strands")

        five = coding_regions.Start.min()
        three = coding_regions.End.max()

        gr_start = gr
        gr_start.End = gr.Start

        gr_end = gr
        gr_end.Start = gr.End

        closest_start = coding_regions.nearest(gr_start)

        # find the starting location on the genome based on coding sequences
        if closest_start.Distance == 0:  # inside a cds
            start_location = gr.Start
        elif closest_start.Start - gr.Start < 0:  # 5' of the closest cds so take the whole cds
            start_location = closest_start.Start
        elif closest_start.Start - gr.Start > 0:  # 3' of the closest cds so need to take the next cds
            start_exon_rank = closest_start.exon_rank + 1
            if start_exon_rank > coding_regions.exon_rank.max():  # start at 3' utr so no protein
                print("gr starts at a 3' utr there is no protein sequence overlap")
            else:
                start_location = coding_regions[coding_regions.exon_rank == start_exon_rank].Start

        # repeat for the end location
        closest_end = coding_regions.nearest(gr_end)

        if closest_end.Distance == 0:  # inside a cds
            end_location = gr.End
        elif closest_end.End - gr.Start > 0:  # 3' of the closest cds so take the whole cds
            end_location = closest_start.End
        elif closest_end.Start - gr.Start < 0:  # 5' of the closest cds so need to take the prevoius cds
            end_exon_rank = closest_end.exon_rank - 1
            if end_exon_rank < coding_regions.exon_rank.min():  # end at 5' utr so no protein
                print("gr ends at a 5' utr there is no protein sequence overlap")
            else:
                end_location = coding_regions[coding_regions.exon_rank == end_exon_rank].End

        # find aa location here strand matters
        if strand == "+":
            start_aa = (start_location - five) / 3
            end_aa = (end_location - five) / 3
        else:
            start_aa = (three - start_location) / 3
            end_aa = (three - end_location) / 3

        columns = [feature_name, feature_name + "_start", feature_name + "_end"]
        domains = self.mart.query(attributes=columns,
                                  filters={"link_ensembl_transcript_stable_id": tx})
        domains = domains.sort_values(by=domains.columns[1])

        starts = domains.iloc[:, 1].to_list()
        ends = domains.iloc[:, 2].to_list()

        intervals = []
        for i in range(len(starts)):
            intervals.append(pd.Interval(starts[i], ends[i], closed="both"))

        affected = pd.Interval(start_aa, end_aa, closed="both")
        overlaps = [interval for interval in intervals if interval.overlaps(affected)]
        start = min([interval.left for interval in overlaps])
        end = max([interval.right for interval in overlaps])
        start_idx = domains.index[domains.iloc[:, 1] == start]
        end_idx = domains.index[domains.iloc[:, 2] == end]

        return domains.iloc[start_idx[0]:end_idx[0], :]
