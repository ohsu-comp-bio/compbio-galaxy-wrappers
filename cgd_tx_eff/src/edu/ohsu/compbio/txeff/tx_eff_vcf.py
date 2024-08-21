'''
Created on Apr 29, 2022

@author: pleyte

Update a VCF with the variant-transcript details.  
'''

import argparse
from collections import defaultdict
from enum import Enum
import logging.config
import sys

import vcfpy

from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig
from edu.ohsu.compbio.txeff.util.tx_eff_csv import TxEffCsv
from edu.ohsu.compbio.txeff.variant import Variant


class TranscriptEffect(Enum):
    '''
    The transcript effects INFO fields that are added to the VCF file.
    
    Most of these INFO fields are ':' delimited parallel arrays where each index in one field is associated
    with the same index in another field.  The TFX_SPLICE and TFX_REFERENCE_CONTEXT fields apply to all 
    transcripts so they only have a single value.  
    ''' 
    TFX_GENE = 'TFX_GENE'
    TFX_TRANSCRIPT = 'TFX_TRANSCRIPT'
    TFX_G_DOT = 'TFX_G_DOT'
    TFX_HGVSC = 'TFX_HGVSC'
    TFX_HGVSP1 = 'TFX_HGVSP1'
    TFX_HGVSP3 = 'TFX_HGVSP3'
    TFX_SPLICE = 'TFX_SPLICE'
    TFX_VARIANT_EFFECT = 'TFX_VARIANT_EFFECT'
    TFX_VARIANT_TYPE = 'TFX_VARIANT_TYPE'
    TFX_BASE_POSITION = 'TFX_BASE_POSITION'
    TFX_PROTEIN_TRANSCRIPT = 'TFX_PROTEIN_TRANSCRIPT'
    TFX_EXON = 'TFX_EXON'
    TFX_AMINO_ACID_POSITION = 'TFX_AMINO_ACID_POSITION'
    TFX_REFERENCE_CONTEXT = 'TFX_REFERENCE_CONTEXT' 

class TxEffVcf(object):
    def __init__(self, tfx_version: str, in_vcf:str, out_vcf:str):
        self.version = tfx_version
        self.logger = logging.getLogger(__name__)
        self.in_vcf_filename = in_vcf
        self.out_vcf_filename = out_vcf
    
    def _read_vcf(self, vcf_filename: str):
        '''
        Read each variant from the VCF file and return a map of Variant[vcf-record]. Each variant must be unique and have only one ALT allele.  
        '''
        vcf_reader = vcfpy.Reader.from_path(vcf_filename)
        variant_dict = {}
        
        for vcf_record in vcf_reader:
            if len(vcf_record.ALT) > 1:
                raise Exception(f'VCF variants must have just one ALT allele: {vcf_record.CHROM}-{vcf_record.POS}-{vcf_record.REF}-{vcf_record.ALT}')
            
            variant = Variant(vcf_record.CHROM, vcf_record.POS, vcf_record.REF, vcf_record.ALT[0].value)
    
            if variant in variant_dict:
                raise Exception(f"Duplicate variant in {vcf_filename}: {variant}")
            else:
                variant_dict[variant] = vcf_record
            
        vcf_reader.close()
        
        return variant_dict, vcf_reader.header

    def _get_transcripts_dict(self, transcripts: list):
        '''
        Return a dictionary where the key is a variant and the value is a list of transcripts 
        ''' 
        transcript_dict = defaultdict(list)
        
        for transcript in transcripts:        
            variant = Variant(transcript.chromosome, transcript.position, transcript.reference, transcript.alt)
            transcript_dict[variant].append(transcript)
            self.logger.debug(f"{variant} has {len(transcript_dict[variant])} transcripts")
    
        return transcript_dict

    def _add_transcripts_to_vcf(self, vcf_variant_dict, transcript_dict):
        '''
        Transform the transcript details into parallel arrays and add them to each VCF row.  
        '''
        vcf_records = []

        for (variant, vcf_record) in vcf_variant_dict.items():
            transcripts = transcript_dict.get(variant)

            if not transcripts:
                self.logger.warning(f"No transcripts found for VCF variant {variant}")
                vcf_records.append(vcf_record)
                continue
            
            # These values are the same for all transcripts so there is only one value. 
            tfx_splice = self._get_any_splicing(transcripts)            
            tfx_reference_context = self._get_and_confirm_single_value('reference_context', lambda x: x.reference_context, transcripts)

            # These values are different for each transcript
            tfx_base_positions = []
            tfx_exons = []
            tfx_genes = []
            tfx_g_dots = []
            tfx_c_dots = []
            tfx_p1 = []
            tfx_p3 = []
            tfx_refseq_transcripts = []
            tfx_variant_effects = []
            tfx_variant_types = []
            tfx_protein_transcripts = []
            tfx_amino_acid_positions = []
            
            for transcript in transcripts:
                tfx_base_positions.append(transcript.hgvs_base_position)
                tfx_exons.append(transcript.exon)
                tfx_genes.append(transcript.hgnc_gene)
                tfx_g_dots.append(transcript.sequence_variant)
                tfx_c_dots.append(transcript.hgvs_c_dot)
                tfx_p1.append(transcript.hgvs_p_dot_one)
                tfx_p3.append(transcript.hgvs_p_dot_three)
                tfx_refseq_transcripts.append(transcript.refseq_transcript)
                tfx_variant_effects.append(transcript.variant_effect)
                tfx_variant_types.append(transcript.variant_type)
                tfx_protein_transcripts.append(transcript.protein_transcript)
                tfx_amino_acid_positions.append(transcript.hgvs_amino_acid_position)
    
            vcf_record.INFO[TranscriptEffect.TFX_SPLICE.value] = [tfx_splice]
            vcf_record.INFO[TranscriptEffect.TFX_REFERENCE_CONTEXT.value] = [tfx_reference_context]
            vcf_record.INFO[TranscriptEffect.TFX_BASE_POSITION.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_base_positions])]
            vcf_record.INFO[TranscriptEffect.TFX_EXON.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_exons])]
            vcf_record.INFO[TranscriptEffect.TFX_GENE.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_genes])]
            vcf_record.INFO[TranscriptEffect.TFX_G_DOT.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_g_dots])]
            vcf_record.INFO[TranscriptEffect.TFX_HGVSC.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_c_dots])]
            vcf_record.INFO[TranscriptEffect.TFX_HGVSP1.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_p1])]
            vcf_record.INFO[TranscriptEffect.TFX_HGVSP3.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_p3])]
            vcf_record.INFO[TranscriptEffect.TFX_AMINO_ACID_POSITION.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_amino_acid_positions])]
            vcf_record.INFO[TranscriptEffect.TFX_TRANSCRIPT.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_refseq_transcripts])]        
            vcf_record.INFO[TranscriptEffect.TFX_VARIANT_TYPE.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_variant_types])]        
            vcf_record.INFO[TranscriptEffect.TFX_PROTEIN_TRANSCRIPT.value] = [':'.join([self._replaceNoneWithEmpty(x) for x in tfx_protein_transcripts])]

            # Replace spaces with underscore because spaces are not allowed in the INFO field.
            vcf_record.INFO[TranscriptEffect.TFX_VARIANT_EFFECT.value] = [':'.join([self._replaceNoneWithEmpty(x).replace(' ', '_') for x in tfx_variant_effects])]
            
            vcf_records.append(vcf_record)
        
        return vcf_records

    def _get_and_confirm_single_value(self, field_name, getter, transcripts: list):
        '''
        Takes a list of trancripts and a getter function for one of the fields in the transcript. This function
        confirms that all the transcripts have the same value for that field, and returns the value. 
        '''
        if not transcripts:
            return None
        
        # Use getter to extract all values and add them to a set 
        values = set(map(getter, transcripts))
        
        if len(values) != 1:
            raise ValueError(f"All values in {field_name} array should be the same but {transcripts[0]} has multiple: {values}")
        
        return values.pop()
    
    def _get_any_splicing(self, trancripts:list):
        '''
        Return any non-empty splicing value from the transcript list. If a variant is splicing then all the transcripts
        identified by Annovar will be splicing. But HGVS/UTA doesn't know about splicing so its trancripts won't have a splicing value.
        This function checks to see if any of the transcripts are splicing and returns the splicing value if one of them is.    
        '''
        for x in trancripts:
            if x.splicing:
                return x.splicing
            
    def _replaceNoneWithEmpty(self, x: str):
        '''
        Return the empty string if x is None, otherwise return x
        '''
        return "" if x is None else str(x)
    
    def _write_vcf(self, vcf_file_name, header: vcfpy.header.Header, vcf_records: list):
        '''
        Take the VCF records that have been updated with transcript effects and write them to file 
        '''
        writer = vcfpy.Writer.from_path(vcf_file_name, header)
        
        for vcf_record in vcf_records:
            writer.write_record(vcf_record)
            
        writer.close()

    def _update_header(self, header):
        '''
        Add the new transcript effect fields to the VCF header 
        '''
        header.add_line(vcfpy.HeaderLine('tfx_commandline', ' '.join(sys.argv)))
        header.add_line(vcfpy.HeaderLine('TfxVersion', self.version))
        
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_BASE_POSITION.value),
                                                ('Number', '.'),
                                                ('Type', 'String'),
                                                ('Description', 'Coding sequence start position.')]))
    
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_EXON.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'Exon number associated with given transcript.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_GENE.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'HGNC gene symbol.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_G_DOT.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'HGVS g-dot nomenclature.')]))    
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_HGVSC.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'HGVS cdot nomenclature.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_HGVSP1.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'HGVS pdot nomenclature, single letter amino acids.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_HGVSP3.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'HGVS pdot nomenclature, three letter amino acids.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_SPLICE.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'Splice site annotation.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_TRANSCRIPT.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'Transcript identifier.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_VARIANT_EFFECT.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'Variant effect annotation.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_VARIANT_TYPE.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'Variant type or location annotation.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_PROTEIN_TRANSCRIPT.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'Protein transcript.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_AMINO_ACID_POSITION.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'Amino acid start position.')]))
        header.add_info_line(vcfpy.OrderedDict([('ID', TranscriptEffect.TFX_REFERENCE_CONTEXT.value),
                                                                ('Number', '.'),
                                                                ('Type', 'String'),
                                                                ('Description', 'Reference context.')]))
        
    
    def create_vcf(self, transcripts: list):
        '''
        Read VCF variants from file, add transcript effects to the INFO fields and write out a new VCF. 
        '''
        # Transform the transcript list into a dictionary 
        transcript_dict = self._get_transcripts_dict(transcripts)
        self.logger.info(f'There are {len(transcript_dict)} distinct variants among {len(transcripts)} transcripts')
        
        # Read all variant rows from the VCF
        vcf_variant_dict, vcf_header = self._read_vcf(self.in_vcf_filename)
        self.logger.info(f'Read {len(vcf_variant_dict)} distinct variants from {self.in_vcf_filename}')
        
        # Combine VCF variants with transcript effects      
        vcf_transcript_records = self._add_transcripts_to_vcf(vcf_variant_dict, transcript_dict)
        assert len(vcf_transcript_records) == len(vcf_variant_dict), f"The number of variants updated with transcript effects is not the same as the number of variants in the VCF: {len(vcf_transcript_records)} != {len(vcf_variant_dict)}"
        
        # Add new INFO fields to VCF header
        self._update_header(vcf_header)
        
        self.logger.info(f"Writing {len(vcf_transcript_records)} variants with transcript effects to VCF file {self.out_vcf_filename}")
        self._write_vcf(self.out_vcf_filename, vcf_header, vcf_transcript_records)
    
def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Produce a VCF containing transcript details generated by tx_eff.hgvs.py')

    parser.add_argument('-i', '--in_vcf', 
                    help='Input VCF', 
                    type=argparse.FileType('r'), 
                    required=True)

    parser.add_argument('-o', '--out_vcf', 
                    help='Output VCF', 
                    type=argparse.FileType('w'), 
                    required=True)

    parser.add_argument('-c', '--in_csv', 
                    help='Input CSV', 
                    type=argparse.FileType('r'), 
                    required=True)
    
    args = parser.parse_args()
    
    return args
    
    args = parser.parse_args()

def _main():
    '''
    main function
    '''
    logging.config.dictConfig(TfxLogConfig().log_config)
    
    args = _parse_args()
    
    logging.info(f"Reading transcripts from {args.in_csv.name}")
    tx_eff_csv = TxEffCsv()    
    transcripts = tx_eff_csv.read_transcripts(args.in_csv.name)
    
    # Combine transcript effects with VCF variants and write to file 
    TxEffVcf(args.in_vcf.name, args.out_vcf.name).create_vcf(transcripts)


if __name__ == '__main__':
    _main()