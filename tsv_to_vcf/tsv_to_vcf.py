#!/usr/bin/env python

### Create a "fake" VCF file based on coordinates inputted by the user.  These are primarily
### used with Sanger results to get quick annotation via SeattleSeq.
### John Letaw 06/01/15

import argparse

class AlleleRestrict(argparse.Action):
    
    def __call__(self, parser, namespace, values, option_string=None):

        bases = set('ACGT')
        if any((b not in bases) for b in list(values.upper())):
            parser.error("Inputs must contain only A, C, G, and T.")

        setattr(namespace, self.dest, values.upper())
        

def writeFakeDoc(chrom, coord):
    return (chrom + ':' + coord + '\t100\t100.00\t100\n')


def main():
    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(dest='chromosome', help='')
    parser.add_argument(dest='position', help='')
    parser.add_argument(dest='ref_allele', action=AlleleRestrict, help='')
    parser.add_argument(dest='alt_allele', action=AlleleRestrict, help='')
    parser.add_argument(dest='gene_name', help='')
    parser.add_argument(dest='output_vcf', help='')
    parser.add_argument(dest='output_gene', help='')
    parser.add_argument(dest='output_doc', help='')
    args = parser.parse_args()

    ### Write the VCF for a single variant.

    handle_vcf = open(args.output_vcf, 'w')
    handle_vcf.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    if ',' in args.chromosome:
        for i in range(len(args.chromosome.split(','))):
            handle_vcf.write(args.chromosome.split(',')[i] + '\t' + args.position.split(',')[i] + '\t.\t' + \
                                 args.ref_allele.split(',')[i] + '\t' + args.alt_allele.split(',')[i] + '\t60\t.\tDB\tGT:AD\t0/1:50,50\n')
    else:
        handle_vcf.write(args.chromosome + '\t' + args.position + '\t.\t' + args.ref_allele + '\t' + args.alt_allele + '\t60\t.\tDB\tGT:AD\t0/1:50,50\n')

    handle_vcf.close()

    ### Write the Gene List based on input from user.

    handle_gene = open(args.output_gene, 'w')
    if ',' in args.gene_name:
        for i in range(len(args.gene_name.split(','))):
            handle_gene.write(args.gene_name(',')[i] + '\n')
    else:
        handle_gene.write(args.gene_name.upper() + '\n')
        handle_gene.write(args.gene_name + '\n')

    handle_gene.close()

    ### Write fake depth of coverage output.

    handle_doc = open(args.output_doc, 'w')
    handle_doc.write("Locus\tTotal_Depth\tAverage_Depth_sample\tDepth_for_normal\n")

    # These are hard-coded in to the writeGenotypes code, so the program will expect to see these in DOC output.
    handle_doc.write(writeFakeDoc("7", "6026775"))
    handle_doc.write(writeFakeDoc("13", "32929387"))
    handle_doc.write(writeFakeDoc("1", "169519049"))
    handle_doc.write(writeFakeDoc("14", "75513883"))
    if ',' in args.chromosome:
        for i in range(len(args.chromosome.split(','))):
            handle_doc.write(writeFakeDoc(args.chromosome.split(',')[i], args.position.split(',')[i]))
    else:
        handle_doc.write(writeFakeDoc(args.chromosome, args.position))
            
    handle_doc.close()


if __name__ == "__main__":
    main()

