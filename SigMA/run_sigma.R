
args <- commandArgs(TRUE)
tumor <- args[1]
platform <- args[2]
lite <- as.logical(as.integer(args[3]))
mva <- as.logical(as.integer(args[4]))
assign <- as.logical(as.integer(args[5]))
output <- args[6]
vcf_dir <- args[7]

print(args)

print("running")

# short example 2 simulated panels for testing the tool
devtools::load_all()
print("loaded")

# data_dir can be replaced by the character containing 
# the directory defined by the user
data_dir <- vcf_dir

print("assigned")

genomes_matrix <- make_matrix(data_dir, file_type = 'vcf')
genomes <- conv_snv_matrix_to_df(genomes_matrix)

print("matrix created")

genome_file = 'vcfs.csv'

write.table(genomes,
            genome_file,
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)

print("genomes writted")

message(paste0('96-dimensional matrix is saved in ', genome_file))
message('Running SigMA')

print(mva)
print(assign)

output_name <- run(genome_file, 
    data = platform, 
    do_assign = assign, 
    do_mva = mva,
    tumor_type = tumor, 
    lite_format = lite)
    
file.rename(output_name, output)
    
print(output)
print("R done")
