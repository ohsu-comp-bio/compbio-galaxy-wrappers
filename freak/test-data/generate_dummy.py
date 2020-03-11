import pandas as pd
import random


# Create references and seed for RNG
bases = set(['T','C','G','A'])
chromosomes = ['chr'+str(c) for c in range(1,22)]+['chrX']
random.seed()

# Generate the dummy reference data
ref_df = pd.DataFrame({'Sample Code': [],
                       'Position Start': [],
                       'Reference Base': [],
                       'Chromosome': [],
                       'Variant Base': []})
for sampleCode in range(100):
    for numVariants in range(random.randint(50,150)):
        idx = random.randint(1,100000000)
        rbase = random.choice([*bases])
        vbase = random.choice([*(bases-set([rbase]))])
        chrm = random.choice(chromosomes)
        row = pd.DataFrame({'Sample Code': [sampleCode+1],
                            'Position Start': [idx],
                            'Reference Base': [rbase],
                            'Chromosome': [chrm],
                            'Variant Base': [vbase]})
        ref_df = ref_df.append(row)
ref_df.loc[:,'Position Start'] = ref_df.loc[:,'Position Start'].astype(int)
ref_df.loc[:,'Sample Code'] = ref_df.loc[:,'Sample Code'].astype(int)

# Sample from the reference data to create dummy freak data
target = random.randint(1,101)
subset_df = ref_df[ref_df['Sample Code']==target].sample(frac=.2).reset_index()
# Alter the chromosome names for the new file
for i,row in subset_df.iterrows():
    print('{}'.format(row['Chromosome'][3:]))
    subset_df.at[i,'Chromosome'] = row['Chromosome'][3:]
print('{}'.format(subset_df.at[1,'Chromosome']))
# Create some false reads
false_df = subset_df.sample(frac=.3)
for i, row in false_df.iterrows():
    used = set([row['Reference Base'],row['Variant Base']])
    false_read = random.choice([*(bases-used)])
    false_df.loc[i,'Variant Base'] = false_read
subset_df = subset_df.append(false_df)

# Rename the columns of the dummy freak output
name_map = {'Position Start': 'start (1-based, inclusive)',
            'Chromosome': 'rname',
            'Variant Base': 'variant',
            'Reference Base': 'reference'}
subset_df = subset_df.rename(columns=name_map)

# Write out the files
ref_df.sample(frac=1).to_csv('data/dummy_ref.tsv',sep='\t',header=True,index=False)
subset_df.sample(frac=1).to_csv(f'data/dummy_freak_{target}.tsv',sep='\t',header=True,index=False)
