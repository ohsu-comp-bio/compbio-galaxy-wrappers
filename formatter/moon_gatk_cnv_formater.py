"""
	Takes a CNV VCF that is formatted with <DEL><DUP> and SNVTYPE=CNV and turns it into a vcf that Moon can read.
"""


import sys


def main(vcf_file, new_vcf_file):
	#open the vcf file
	vcf = open(vcf_file, 'r')
	#open the new file to write to
	new_vcf = open(new_vcf_file, 'w')
	#loop over each of the records in the vcf
	for line in vcf:
		#if its a header line just write it back to the new vcf
		if line.startswith('#'):
			new_vcf.write(line)
		else:
			#Split the line into workable entries
			line_array = line.split('\t')
			#Get the sample entry
			sample_entry = line_array[9]
			#Get the type of the entry from the first character of the sample entry
			if sample_entry.startswith('1'):
				CNV_type = 'DEL'
			elif sample_entry.startswith('2'):
				CNV_type = 'DUP'
			else:
				CNV_type = 'NONE'
			#Change the ALT
			line_array[4] = '<' + CNV_type + '>'
			
			#get the info entry
			info = line_array[7]
			#get the svtype entry 
			info_array = info.split(';')
            
			#change the svtype entry
			if (len(info_array) > 6):
				info_array[6] = 'SVTYPE=' + CNV_type
			else:
				sv_type = 'SVTYPE=' + CNV_type
				info_array.append(sv_type)
			#change the info entry
			corrected_info = ''
			first = True
			for entry in info_array:
				if first:
					corrected_info = entry
					first = False
				else:
					corrected_info += ';' + entry
			#add the corrected info back into the line_array
			line_array[7] = corrected_info
			
			if not CNV_type == "NONE":
				#write the fixed line array to the new vcf
				first = True
				for entry in line_array:
					if first:
						new_vcf.write(entry)
						first = False
					else:
						new_vcf.write('\t' + entry)
	vcf.close()
	new_vcf.close()


if __name__ == "__main__":
    args = sys.argv
    main(args[1], args[2])