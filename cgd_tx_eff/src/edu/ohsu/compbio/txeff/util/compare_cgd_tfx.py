'''
Created on Jul 1, 2022

Takes as input a VCF that has been updated with transcript effects, and a variant export CSV, and 
joins the two so we can compare the nomenclature used in transcript effects with what is currently in CGD. 

@author: pleyte
'''
import argparse
import logging
import csv
import vcfpy
from edu.ohsu.compbio.txeff.tx_eff_vcf import TranscriptEffect

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

stream_handler = logging.StreamHandler()
logging_format = '%(levelname)s: [%(filename)s:%(lineno)s - %(funcName)s()]: %(message)s'

stream_format = logging.Formatter(logging_format)
stream_handler.setFormatter(stream_format)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)

def find_in_vcf(chromosome: str, position_start: int, position_end: int, ref: str, alt: str, vcf_reader: vcfpy.reader.Reader):
    '''
    Return the variant context matching the given genotype
    '''
    chromosome_number = chromosome[3:] if chromosome.startswith("chr") else chromosome
    if position_start > position_end:
        logger.warning("Unable to search for variant because {position_start}>{position_end}: {chromosome}-{position_start}-{ref}-{alt}")
        return None

    vcf_records = [] 
    
    search_string = f"{chromosome_number}:{position_start}-{position_end}"
    
    for record in vcf_reader.fetch(search_string):
        logger.debug(f"Comparing {chromosome_number}-{position_start}-{ref}-{alt} with {record.CHROM}-{record.POS}-{record.REF}-{record.ALT[0].value}")
        if position_start == record.POS and ref == record.REF and alt == record.ALT[0].value:
            vcf_records.append(record)
    
    if(len(vcf_records) == 0):
        logger.debug(f"The VCF does not have a row matching {chromosome_number}-{position_start}-{ref}-{alt}")
    elif(len(vcf_records) > 1):
        logger.error(f"The VCF has more than one row matching {chromosome_number}-{position_start}-{ref}-{alt}. Skipping.")
    elif(len(record.ALT) > 1):
        logger.error(f"The matching VCF record has more than one ALT value: {chromosome_number}-{position_start}-{ref}-{alt} ALT={record.ALT}")
    else:
        logger.debug(f"Found {chromosome_number}-{position_start}-{ref}-{alt}")
        return vcf_records[0]

    return None

def process(vcf_filename: str, csv_export_filename: str, csv_out_filename: str):
    out_fields = ['genomic_variant_id', 'chromosome', 'position', 'reference', 'alt',
                  'cgd_variant_effect', 'cgd_variant_type', 'cgd_hgvs_amino_acid_position', 'cgd_hgvs_base_position', 
                  'cgd_exon', 'cgd_hgnc_gene', 'cgd_primary_gene', 'cgd_hgvs_c_dot', 'cgd_hgvs_p_dot_one', 'cgd_hgvs_p_dot_three', 
                  'cgd_splicing', 'cgd_refseq_transcript', 'cgd_protein_transcript',
                  'tfx_variant_effect', 'tfx_variant_type', 'tfx_hgvs_amino_acid_position', 'tfx_hgvs_base_position', 
                  'tfx_exon', 'tfx_hgnc_gene', 'tfx_hgvs_c_dot', 'tfx_hgvs_p_dot_one', 'tfx_hgvs_p_dot_three', 
                  'tfx_splicing', 'tfx_refseq_transcript', 'tfx_protein_transcript']
    
    logger.info(f"Reading {vcf_filename}")
    vcf_reader = vcfpy.Reader.from_path(vcf_filename)
    
    logger.info(f"Reading {csv_export_filename}")
    
    with open(csv_export_filename) as csv_export_file, open(csv_out_filename, "w") as csv_out_file:
        # Open the variant export csv 
        export_reader = csv.DictReader(csv_export_file)

        # Open the output csv
        csv_out_writer = csv.writer(csv_out_file)
        csv_out_writer.writerow(out_fields)
        
        # Iterate over each entry in the export file 
        for export_row in export_reader:
            logger.debug(f"Looking for {export_row['chromosome']}-{export_row['position_start']}-{export_row['reference_base']}-{export_row['variant_base']}")
            
            vc = find_in_vcf(export_row['chromosome'], int(export_row['position_start']), int(export_row['position_end']), export_row['reference_base'], export_row['variant_base'], vcf_reader)
            if (vc == None):
                logger.info(f"Unable to find genotype {export_row['chromosome']}-{export_row['position_start']}-{export_row['reference_base']}-{export_row['variant_base']} in vcf")
                continue
            
            if vc.INFO.get(TranscriptEffect.TFX_VARIANT_EFFECT.value) == None:
                logger.info(f"Variant does not have variant effects: {export_row['chromosome']}-{export_row['position_start']}-{export_row['reference_base']}-{export_row['variant_base']}")
                continue
            
            # Write row(s) related to this genotype to the csv 
            write(csv_out_writer, export_row, vc)

    logger.info(f"Wrote {csv_out_filename}")

def write(csv_writer, export_row, vc: vcfpy.record.Record):
    '''
    '''    
    csv_writer.writerow([export_row['genomic_variant_id'],
                         export_row['chromosome'], export_row['position_start'], export_row['reference_base'], export_row['variant_base'],
                         export_row['protein_variant_type'],                          # 'cgd_variant_effect' - annovar variant_effect is cgd's protein variant type
                         export_row['genomic_variant_type'],                          # 'cgd_variant_type'
                         export_row['amino_acid_position'],                           # 'cgd_hgvs_amino_acid_position' 
                         export_row['base_pair_position'],                            # 'cgd_hgvs_base_position', 
                         export_row['exon'],                                          # 'cgd_exon'
                         export_row['genes'],                                         # 'cgd_hgnc_gene'
                         export_row['primary_gene'],                                   # 'cgd_primary_gene'
                         export_row['genotype_cdna'],                                 # 'cgd_hgvs_c_dot'
                         export_row['genotype_amino_acid_onel'],                      # 'cgd_hgvs_p_dot_one' 
                         export_row['genotype_amino_acid_threel'],                    # 'cgd_hgvs_p_dot_three' 
                         export_row['splice_site_relationship'],                      # 'cgd_splicing'
                         export_row['cdna_transcript'],                               # 'cgd_refseq_transcript'
                         export_row['protein_transcript'],                           # 'cgd_protein_transcript'
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_VARIANT_EFFECT.value]]),         # 'tfx_variant_effect'
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_VARIANT_TYPE.value]]),           # 'tfx_variant_type', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_AMINO_ACID_POSITION.value]]),    # 'tfx_hgvs_amino_acid_position', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_BASE_POSITION.value]]),          # 'tfx_hgvs_base_position', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_EXON.value]]),                   # 'tfx_exon', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_GENE.value]]),                   # 'tfx_hgnc_gene', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_HGVSC.value]]),                  # 'tfx_hgvs_c_dot', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_HGVSP1.value]]),                 # 'tfx_hgvs_p_dot_one', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_HGVSP3.value]]),                 # 'tfx_hgvs_p_dot_three', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_SPLICE.value]]),                 # 'tfx_splicing', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_TRANSCRIPT.value]]),             # 'tfx_refseq_transcript', 
                        ' / '.join([x for x in vc.INFO[TranscriptEffect.TFX_PROTEIN_TRANSCRIPT.value]])])    # 'tfx_protein_transcript'
    
    

def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Join a VCF containing transcript effects with a CSV from variant export')

    parser.add_argument('-v', '--vcf',  
                        help='Input VCF.bgz containing transcript effects (index must also exist)',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('-c', '--export_csv',  
                        help='Variant export CSV (use variant_transcripts.sql to generate the csv)',
                        type=argparse.FileType('r'),
                        required=True)
        
    parser.add_argument('-o', '--out_file',  
                        help='CSV file to create',
                        type=argparse.FileType('w'),
                        required=True)    
    
    args = parser.parse_args()
    
    return args

def _main():
    '''
    main function
    '''
    args = _parse_args()
    process(args.vcf.name, args.export_csv.name, args.out_file.name)

if __name__ == '__main__':
    _main()