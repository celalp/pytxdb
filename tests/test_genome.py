from Bio.Seq import Seq
import pytest
from pytxdb import genome
from pyranges import *
import sys
from sqlalchemy import create_engine
import pybiomart as biomart

db=create_engine("sqlite:////hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/test.db")
mart=biomart.Dataset(name='hsapiens_gene_ensembl',
                    host='grch37.ensembl.org')
fasta="/hpf/largeprojects/ccmbio/yliang/test_place/pytxdb_test/selected_chr22.fa"

mygenome=genome.Genome(db, fasta, mart)

class TestGenomeclass:
    
    def test_gene_coordinate_search(self):
        target_generegion = mygenome.genes(names=['ENSG00000128272'])
        assert target_generegion.Start[0] == 39915699
        assert target_generegion.End[0] == 39918691
        assert target_generegion.Strand[0] == '+'

    def test_transcript_search(self):
        target_transcript = mygenome.transcripts(gene_names=['ENSG00000128272'])
        assert 'ENST00000404241' in target_transcript.id.values
        assert 'ENST00000337304' in target_transcript.id.values
        assert 'ENST00000396680' in target_transcript.id.values
    
    def test_exons_search(self):
        #is the first exon for this transcript equal to ENSE00001552405
        target_exons = mygenome.exons(names=['ENST00000404241'])
        assert 'ENSE00001552405' in target_exons.exon_name[target_exons.rank==1].values
        assert 'ENSE00001558792' in target_exons.exon_name[target_exons.rank==2].values
        assert 'ENSE00001525916' in target_exons.exon_name[target_exons.rank==3].values
        assert 'ENSE00001561335' in target_exons.exon_name[target_exons.rank==4].values
    
    def test_cds_search(self):
        target_cds = mygenome.cds(names=['ENST00000404241'], level="transcript")
        assert target_cds.Start.tolist() == [39917450, 39917777]
    
    def test_intron_search(self):
        target_intron = mygenome.introns(names=['ENST00000404241'], level="transcript")
        assert target_intron.Start.tolist() == [39915807, 39916757, 39917677]
        
    def test_3utr_search(self):
        target_3utr = mygenome.three_utr(names=['ENST00000404241'], level="transcript")
        assert target_3utr.Start[0] == 39918607
        assert target_3utr.End[0] == 39918688
        
    def test_5utr_search(self):
        target_5utr = mygenome.five_utr(names=['ENST00000337304'], level="transcript")
        assert target_5utr.Start[0] == 39916568
        assert target_5utr.End[0] == 39917450
    
    def test_domain_search(self):
        gr = PyRanges(chromosomes='22', starts=[36677327], ends=[36784063],strands=("-"))
        domain_df = mygenome.domains(gr=gr,tx='ENST00000216181',feature_name='pfam')
        assert 'PF00063' in domain_df[0]['Pfam ID'].values
        assert 'PF00612' in domain_df[0]['Pfam ID'].values
        assert 'PF01576' in domain_df[0]['Pfam ID'].values
        assert 'PF02736' in domain_df[0]['Pfam ID'].values

    def test_get_sequence(self):
        gr = PyRanges(chromosomes="22", starts=[39916568], ends=[39916578],strands=("+"))
        target_sequence = mygenome.get_sequence(gr=gr, type="nuc", concat=False)[0]
    
        assert target_sequence == Seq('TTTCTACTTT')

    def test_gene_search(self):
        target_gene = mygenome.search('^MYH9$|^ADM2$',where='hgnc_symbol', regex=True)
        assert target_gene == ['ENSG00000100345', 'ENSG00000128165']
    
