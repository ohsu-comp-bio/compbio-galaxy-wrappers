import argparse

VERSION = '0.0.5'

def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('gff', help='')
    parser.add_argument('vcf', help='')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('gene_list', help='')
    parser.add_argument('--ordered_test', '-o')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def filter_genes(genelist_filepath, ot_input):

    gene_list = open(genelist_filepath, 'r')
    gene_filter = []

    match = False
    for line in gene_list:
        line_array = line.split('\t')
        ordered_test = line_array[1]

        if ot_input == ordered_test:
            match = True
            gene_filter.append(line_array[2])

    if match == False:
        raise NameError('Please enter proper ordered test name.')

    return gene_filter


def annotate_vcf(gff_filepath, vcf_filepath, new_filepath, genes=[]):
    gff = open(gff_filepath, 'r')
    vcf = open(vcf_filepath, 'r')
    new_vcf = open(new_filepath, 'w')

    # Init dictionary of GENE:[chr, start pos, end pos] from GFF
    gene_dict = {}

    # loop over each of the records in the vcf
    for line in gff:
        if line.startswith('#'):
            continue
        else:
            # Split the line into workable entries
            line_array = line.split('\t')
            if line_array[2] == 'region':
                info = line_array[8].split(';')[3]
                current_chr = info.partition('chromosome=')[2]
            elif line_array[2] == 'gene':
                start = line_array[3]
                end = line_array[4]
                gene = line_array[8].split(';')[1].partition('Name=')[2]
                gene_dict[gene] = [current_chr, start, end]

    # Insert ##INFO description for Gene annotation
    header_insert = False
    for line in vcf:
        # If it's a header line, just write it to new VCF
        if line.startswith('#'):
            new_vcf.write(line)
            if line.startswith('##INFO') and header_insert is False:
                new_vcf.write(
                    '##INFO=<ID=Gene,Number=.,Type=String,Description="Genes corresponding to coordinates of this variant.">\n')
                header_insert = True
        else:
            # Get INFO entry
            line_array2 = line.split('\t')
            info_array = line_array2[7].split(';')
            info_array_check = len(info_array)

            # iterate through .gff genes
            for k, v in gene_dict.items():
                if v[0] == line_array2[0] and int(line_array2[1]) >= int(v[1]) and int(line_array2[1]) <= int(v[2]) and len(info_array) == info_array_check:
                    info_array.append('Gene=' + k)
                elif v[0] == line_array2[0] and int(line_array2[1]) >= int(v[1]) and int(line_array2[1]) <= int(v[2]) and len(info_array) > info_array_check:
                    info_array[-1] = str(info_array[-1]) + ',' + k

            # Filter if not in gene list
            gene_in_list = False
            if info_array[-1][0:4] == 'Gene' and len(genes) != 0:
                for g in info_array[-1].split(','):
                    if g[0:5] == 'Gene=' and g[5:] in genes:
                        gene_in_list = True
                    elif g in genes:
                        gene_in_list = True

            if gene_in_list is False and len(genes) != 0:
                continue

            # Filter out non-annotated lines by skipping it if gene info isn't available
            if info_array[-1][0:4] != 'Gene' and len(genes) != 0:
                continue

            first = True
            for entry in info_array:
                if first:
                    corrected_info = entry
                    first = False
                else:
                    corrected_info += ';' + entry
            line_array2[7] = corrected_info

            # write the updated line array to the new vcf
            first = True
            for entry in line_array2:
                if first:
                    new_vcf.write(entry)
                    first = False
                else:
                    new_vcf.write('\t' + entry)

    vcf.close()
    new_vcf.close()
    return gene_dict.keys()

def main():

    args = supply_args()

    if args.ordered_test != None:
        genes = filter_genes(args.gene_list, args.ordered_test)
        gff_genes = list(annotate_vcf(args.gff, args.vcf, args.outfile, genes))
        print('Ordered test genes: ', len(genes))
        print('GFF genes: ', len(gff_genes))
        print('Number of common genes: ', len(set(genes) & set(gff_genes)))
        print('Common: ', set(genes) & set(gff_genes))
        # Stop tool and print missing genes, if applicable
        uncommon = []
        for i in set(genes):
            if i not in set(gff_genes):
                uncommon.append(i)
        if len(uncommon) > 0:
            raise AttributeError(uncommon, ' is/are missing from reference.')
    else:
        annotate_vcf(args.gff, args.vcf, args.outfile)

if __name__ == '__main__':
    main()