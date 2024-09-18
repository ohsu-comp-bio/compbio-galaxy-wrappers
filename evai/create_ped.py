# Creates trio pedigree file for eVai

import argparse

VERSION = '0.0.1'

parser = argparse.ArgumentParser(description='')
parser.add_argument('--output_ped','-o')
parser.add_argument('family_id', help='')
parser.add_argument('proband_id', help='')
parser.add_argument('proband_sex', help='')
parser.add_argument('father_id', help='')
parser.add_argument('mother_id', help='')
parser.add_argument('--father_ph', '-f')
parser.add_argument('--mother_ph', '-m')

args = parser.parse_args()

# Optional inputs for parental phenotypes
if args.father_ph is None or args.father_ph == '':
    fp = '1'
else:
    fp = args.father_ph
if args.mother_ph is None or args.mother_ph == '':
    mp = '1'
else:
    mp = args.mother_ph

proband = [args.family_id, args.proband_id, args.father_id, args.mother_id, args.proband_sex, '2']
father = [args.family_id, args.father_id, '0', '0', '1', fp]
mother = [args.family_id, args.mother_id, '0', '0', '2', mp]

with open(args.output_ped, 'w') as f:
    f.write('\t'.join(proband) + '\n')
    f.write('\t'.join(father) + '\n')
    f.write('\t'.join(mother))