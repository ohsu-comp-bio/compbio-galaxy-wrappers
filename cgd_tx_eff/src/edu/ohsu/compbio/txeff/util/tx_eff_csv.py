'''
Created on Aug. 24, 2022

A class that can read and write transcripts to/from CSV

@author: pleyte
'''
import csv
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript

class TxEffCsv(object):
    '''
    Read and write VariantTranscript objects to/from CSV
    '''

    def read_transcripts(self, in_file):
        '''
        Read the transcripts that have been written to a CSV file  
        '''
        transcripts = []
    
        with open(in_file) as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                transcript = VariantTranscript(row['chromosome'], row['position'], row['reference'], row['alt'])
                transcript.hgvs_amino_acid_position = row['hgvs_amino_acid_position']
                transcript.variant_effect = row['variant_effect']
                transcript.variant_type = row['variant_type']
                transcript.hgvs_amino_acid_position = self._noneIfEmpty(row['hgvs_amino_acid_position'])
                transcript.hgvs_base_position = self._noneIfEmpty(row['hgvs_base_position'])
                transcript.exon = row['exon']
                transcript.hgnc_gene = row['hgnc_gene']
                transcript.hgvs_c_dot = row['hgvs_c_dot']
                transcript.hgvs_p_dot_one = row['hgvs_p_dot_one']
                transcript.hgvs_p_dot_three = row['hgvs_p_dot_three']
                transcript.splicing = row['splicing']
                transcript.refseq_transcript = row['refseq_transcript']
                transcript.protein_transcript = row['protein_transcript']
                
                transcripts.append(transcript)
       
        return transcripts

    def write_transcripts(self, out_file, transcripts: list):
        '''
        Write a list of VariantTranscript records to csv file. 
        '''
        fields = ['chromosome', 'position', 'reference', 'alt', 
          'variant_effect', 'variant_type', 'hgvs_amino_acid_position', 'hgvs_base_position', 
          'exon', 'hgnc_gene', 'hgvs_c_dot', 'hgvs_p_dot_one', 'hgvs_p_dot_three', 
          'splicing', 'refseq_transcript', 'protein_transcript']
    
        with open(out_file, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(fields)
            
            for rec in transcripts:                     
                if(type(rec) == VariantTranscript.__class__):
                    protein_transcript = rec.protein_transcript
                else: 
                    protein_transcript  = ''
                               
                csv_writer.writerow([rec.chromosome, rec.position, rec.reference, rec.alt, 
                                    rec.variant_effect, rec.variant_type, rec.hgvs_amino_acid_position, rec.hgvs_base_position,
                                    rec.exon, rec.hgnc_gene, rec.hgvs_c_dot, rec.hgvs_p_dot_one, rec.hgvs_p_dot_three,
                                    rec.splicing, rec.refseq_transcript, protein_transcript])

    def _noneIfEmpty(self, x):
        '''
        Return the empty string if x is None, otherwise return x
        '''
        return "" if x is None else str(x)