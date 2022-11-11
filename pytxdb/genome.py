import os
import numpy as np
import pandas as pd
import pysam
import sqlalchemy as sql
from Bio.Seq import Seq
from sqlalchemy.orm import Session
from pytxdb.utils import check_results, GRanges

class Genome:
    def __init__(self, db, fasta=None, mart=None):
        """
        get the database connection, and fasta connection for getting sequences
        :param db: genome database, this is a connection
        :param fasta: fasta file for getting sequences optional.
        :param mart pybiomart connection optional
        """
        self.db = db
        if fasta is not None:
            if not os.path.isfile(fasta):
                raise FileNotFoundError("{} does not exist".format(fasta))
            else:
                self.fasta = pysam.Fastafile(fasta)
        self.mart = mart
        self.session = Session(self.db)
        self.metadata = sql.MetaData(self.db)
        self.metadata.reflect(bind=self.db)

    def search(self, term=None, where=None, regex=False, return_fields=False, gene=True):
        """
        searches for ensembl id based on some information like gene name can use regexes
        :param term: what to search for
        :param where: optional search field
        :param regex: is this a regex default False
        :param gene: search gene annotations table othewise transcript annotations
        :return: returns a list of strings (ens ids)
        """
        if gene:
            annot_table = sql.Table("gene_annotations", self.metadata, autoload=True, autoload_with=self.db)
            id="ensembl_gene_id"
        else:
            annot_table = sql.Table("transcript_annotations", self.metadata, autoload=True, autoload_with=self.db)
            id = "ensembl_transcript_id"

        if return_fields:
            print("Available fields are: {}".format(",".join(annot_table.columns.keys())))
        else:

            query=sql.select(annot_table)
            annots_table=self.session.execute(query).fetchall()
            annots=pd.DataFrame(annots_table)

            annots.columns=[desc["name"] for desc in query.column_descriptions]

            mask_array = np.zeros_like(annots, dtype=bool)
            if where is None:
                for i in range(len(annots.columns)):
                    mask_array[:, i] = annots.iloc[:, i]. \
                        str.contains(term, case=False, regex=regex, na=False)
                mask = np.sum(mask_array, axis=1) > 0

            else:
                mask = annots[where].str.contains(term, case=False, regex=regex, na=False)

            results = annots[id][mask].drop_duplicates()
            return results.to_list()

    def genes(self, names=None, df=False, add_annotations=False):
        """
        get gene coordinates
        :param names: list of ens ids if none everything, this needs to be an iterable even
        if only one thing is requested
        :param df: return a dataframe or pyranges
        :param additional_features: add additional features to the search results
        :return: returns a dataframe or pyranges object
        """
        table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)
        if names is None:
            query = sql.select(table)
        else:
            query = sql.select(table).filter(table.c.id.in_(names))

        if add_annotations:
            annot_table=sql.Table("gene_annotations", self.metadata,
                                  autoload=True, autoload_with=self.db)
            query=query.join(annot_table, table.c.id == annot_table.c.ensembl_gene_id).\
                add_columns(*[annot_table.c[col] for col in annot_table.columns.keys()])

        results = check_results(self.session.execute(query).fetchall())

        columns = [desc["name"] for desc in query.column_descriptions]

        dat = pd.DataFrame(results)
        dat.columns = columns

        if df:
            return dat
        else:
            gr = GRanges(chromosomes=dat.iloc[:, 1],
                                            starts=dat.iloc[:, 2],
                                            ends=dat.iloc[:, 3],
                                            strands=dat.iloc[:, 4])
            gr.id = dat["id"]
            return gr

    def transcripts(self, names=None, gene_names=None, df=False, add_annotations=False):
        """
        return transcripts
        :param names: list of ens ids if none everything, this needs to be an iterable even if
        there is only one transcript requested
        :param df: return a df otherwise pyranges
        :param add_annotations: join with transcript annotations table
        :return: a df or pyranges
        """
        tx_table = sql.Table("transcripts", self.metadata, autoload=True, autoload_with=self.db)
        gene_table = sql.Table("genes", self.metadata, autoload=True, autoload_with=self.db)
        # need to get chrom and strand
        subq = sql.select(gene_table.c.id, gene_table.c.chrom, gene_table.c.strand).subquery()

        query = sql.select(tx_table, subq.c.chrom, subq.c.strand). \
            join(subq, tx_table.c.gene == subq.c.id)

        if names is not None:
            query = query.filter(tx_table.c.id.in_(names))

        if gene_names is not None:
            query = query.filter(tx_table.c.gene.in_(gene_names))

        if add_annotations:
            annot_table=sql.Table("transcript_annotations", self.metadata,
                                  autoload=True, autoload_with=self.db)
            query=query.join(annot_table, tx_table.c.id == annot_table.c.ensembl_transcript_id).\
                add_columns(*[annot_table.c[col] for col in annot_table.columns.keys()])

        results = check_results(self.session.execute(query).fetchall())
        dat = pd.DataFrame(results)
        dat.columns = [desc["name"] for desc in query.column_descriptions]

        if df:
            return dat
        else:
            gr = GRanges(chromosomes=dat.iloc[:, 5],
                                            starts=dat.iloc[:, 2],
                                            ends=dat.iloc[:, 3],
                                            strands=dat.iloc[:, 6])
            gr.tx_id = dat.id
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
        dat.columns =[desc["name"] for desc in query.column_descriptions]
        dat = dat.drop(columns="id")

        if df:
            return dat
        else:
            gr = GRanges(chromosomes=dat.iloc[:, 6],
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
        dat.columns = [desc["name"] for desc in query.column_descriptions]
        dat = dat.drop(columns=["id"])

        if df:
            return dat
        else:
            gr = GRanges(chromosomes=dat.iloc[:, 4],
                                            starts=dat.iloc[:, 1],
                                            ends=dat.iloc[:, 2],
                                            strands=dat.iloc[:, 5])
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
        dat.columns = [desc["name"] for desc in query.column_descriptions]
        dat = dat.drop(columns="id")

        if df:
            return dat
        else:
            gr = GRanges(chromosomes=dat.iloc[:, 4],
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
        dat.columns = [desc["name"] for desc in query.column_descriptions]

        if df:
            return dat
        else:
            gr = GRanges(chromosomes=dat.iloc[:, 5],
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
        dat.columns = [desc["name"] for desc in query.column_descriptions]

        if df:
            return dat
        else:
            gr = GRanges(chromosomes=dat.iloc[:, 5],
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

    def domains(self, gr, feature_name, tx, quiet=True):
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

        # pranges returns strand specific
        coding_regions = self.cds(names=[tx], level="transcript", df=False)
        coding_regions.exon_rank = [int(rank) for rank in coding_regions.exon_rank]
        coding_regions.len = coding_regions.End - coding_regions.Start

        strand = coding_regions.Strand.unique()[0]

        if strand != gr.Strand[0]:
            raise ValueError("region provided and the transcript are on different strands")

        coding_regions.relative_end = np.cumsum(coding_regions.len)
        coding_regions.relative_start = pd.concat([pd.Series([1]), coding_regions.relative_end[:-1] + 1])

        gr_start = gr.copy()
        gr_start.End = gr.Start

        gr_end = gr.copy()
        gr_end.Start = gr.End

        closest_start = gr_start.nearest(coding_regions)
        #intialize variables here and they will be overwritten later
        message=None
        start_aa=None
        end_aa=None
        domains=None

        if int(closest_start.End_b) - int(gr.Start) <= 0:  # 3' of the closest cds need to get the next one
            if strand == "+":
                start_exon_rank = int(closest_start.exon_rank) + 1  # get the next one 3' of the exon
                if start_exon_rank > coding_regions.exon_rank.max():  # start at 3' utr so no protein
                    message="gr starts at a 3' utr there is no protein sequence overlap"
                    domains = None
                    start_aa = 0
                    end_aa = 0
                else:
                    start_location = int(coding_regions[coding_regions.exon_rank == start_exon_rank].Start)
            else:
                start_exon_rank = int(closest_start.exon_rank) - 1  # because exon rank decreases as we go 3'
                if start_exon_rank < coding_regions.exon_rank.min():  # start at 5' utr so no protein
                    message="gr ends at a 5' utr there is no protein sequence overlap"
                    domains = None
                    start_aa = 0
                    end_aa = 0
                else:
                    start_location = int(coding_regions[coding_regions.exon_rank == start_exon_rank].Start)
        elif int(closest_start.Start_b) - int(gr.Start) >= 0:  # 5' of the closest cds get the whole thing
            start_location = int(closest_start.Start_b)
            start_exon_rank = int(closest_start.exon_rank)
        else:  # inside the cds
            start_exon_rank = int(closest_start.exon_rank)
            start_location = int(gr.Start)  # because inside the exon

        # repeat for the end location
        closest_end = gr_end.nearest(coding_regions)

        if int(closest_end.Start_b) - int(gr.End) >= 0:  # skipping a cds and the closest one is the next cds
            if strand == "+":
                end_exon_rank = int(closest_end.exon_rank) - 1
                if end_exon_rank < coding_regions.exon_rank.min():  # end at 5' utr so no protein
                    message="gr ends at a 5' utr there is no protein sequence overlap"
                    domains = None
                    start_aa = 0
                    end_aa = 0
                else:
                    end_location = int(coding_regions[coding_regions.exon_rank == end_exon_rank].End)
            else:
                end_exon_rank = int(closest_end.exon_rank) + 1
                if end_exon_rank > coding_regions.exon_rank.max():  # end at 5' utr so no protein
                    message="gr ends at a 5' utr there is no protein sequence overlap"
                    domains = None
                    start_aa = 0
                    end_aa = 0
                else:
                    end_location = int(coding_regions[coding_regions.exon_rank == end_exon_rank].End)
        elif int(closest_end.End_b) - int(gr.End) <= 0:  # closest cds is the one we are skipping
            end_exon_rank = int(closest_end.exon_rank)
            end_location = int(closest_end.End_b)
        else:  # inside the cds
            end_exon_rank = int(closest_end.exon_rank)
            end_location = int(gr.End)

        if start_aa==0 and end_aa==0 and domains is None:
            if not quiet:
                print(message)
            return domains, start_aa, end_aa

        # find aa location here strand matters
        start = coding_regions[coding_regions.exon_rank == start_exon_rank]
        end = coding_regions[coding_regions.exon_rank == end_exon_rank]
        if strand == "+":
            relative_start = start_location - start.Start + start.relative_start
            relative_end = end_location - end.Start + end.relative_start
            start_aa = float(relative_start / 3)
            end_aa = float(relative_end / 3)
        else:
            relative_end = start.End - start_location + start.relative_start
            relative_start = end.End - end_location + end.relative_start
            start_aa = float(relative_start / 3)
            end_aa = float(relative_end / 3)

        # edge case, if they are in the same intron
        if (start_aa > end_aa) and (start_exon_rank == end_exon_rank):
            message="No affected domains found for {}".format(tx)
            domains=None
        else:
            columns = [feature_name, feature_name + "_start", feature_name + "_end"]
            domains = self.mart.query(attributes=columns,
                                      filters={"link_ensembl_transcript_stable_id": tx})
            if pd.isnull(domains.iloc[0, 0]):
                message="No domains found with {} in {}".format(feature_name, tx)
                domains=None
            else:
                domains = domains.sort_values(by=domains.columns[1])
                starts = domains.iloc[:, 1].to_list()
                ends = domains.iloc[:, 2].to_list()

                intervals = []
                for i in range(len(starts)):
                    intervals.append(pd.Interval(starts[i], ends[i], closed="both"))

                affected = pd.Interval(start_aa, end_aa, closed="both")
                overlaps = [interval for interval in intervals if interval.overlaps(affected)]
                if len(overlaps) > 0:
                    start = min([interval.left for interval in overlaps])
                    end = max([interval.right for interval in overlaps])
                    start_idx = domains.index[domains.iloc[:, 1] == start]
                    end_idx = domains.index[domains.iloc[:, 2] == end] + 1
                    domains= domains.iloc[start_idx[0]:end_idx[0], :]
                else:
                    message="No affected domains found for {}".format(tx)
                    domains=None

        if not quiet and message is not None:
            print(message)
        return domains, start_aa, end_aa