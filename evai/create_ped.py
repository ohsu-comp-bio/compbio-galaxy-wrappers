# Creates trio pedigree file for eVai

import argparse

VERSION = '0.0.1'

parser = argparse.ArgumentParser(description='')
parser.add_argument('family_id', help='')
parser.add_argument('proband_id', help='')
parser.add_argument('proband_sex', help='')
parser.add_argument('proband_ph', help='')
parser.add_argument('father_id', help='')
parser.add_argument('father_ph', help='')
parser.add_argument('mother_id', help='')
parser.add_argument('mother_ph', help='')

args = parser.parse_args()

proband = [args.family_id, args.proband_id, args.father_id, args.mother_id, args.proband_sex, args.proband_ph]
father = [args.family_id, args.father_id, '0', '0', '1', args.father_ph]
mother = [args.family_id, args.mother_id, '0', '0', '2', args.mother_ph]

with open('evai.ped', 'w') as f:
    f.write('\t'.join(proband) + '\n')
    f.write('\t'.join(father) + '\n')
    f.write('\t'.join(mother))