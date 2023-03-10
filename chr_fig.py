import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import io
import sys
from tqdm import tqdm

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not '##' in l]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        header=0,
        sep='\t'
    )

vcf_address=sys.argv[1]
out_name=sys.argv[2]

df = read_vcf(vcf_address)
length_ls = pd.read_table('./stacks/catalog.chrs.tsv')
df = df.loc[:, ['#CHROM','POS']]

tmp = list(df['#CHROM'])
chromosome = []
for a in tmp:
    if not a in chromosome:
        chromosome.append(a)

fig = plt.figure(figsize=(30, 4.5 * len(chromosome)))
count = 0
pbar = tqdm(total=10000)
n=int(10000/len(chromosome))
for chrom in chromosome:
    pbar.set_description("Processing %s:" % chrom)
    #     fig, ax1 = plt.subplots(1, figsize=(20, 3))
    # chrom=chromosome[1]
    # print('Dealing with chromosome: {}'.format(chrom))
    count += 1
    ax = fig.add_subplot(20, 3, count)

    data = df[df['#CHROM'] == chrom]
    ls = list(data['POS'])
    length = int(length_ls[length_ls['# Chrom'] == chrom]['Length'])

    features = []

    for ind in data.index:
        pos = int(data.loc[ind, 'POS'])
        features.append(GraphicFeature(start=pos, end=pos, strand=+1, color="#ffd700",
                                       label=str(pos)))
    record = GraphicRecord(sequence_length=length, features=features, ticks_resolution=length / 4)
    ax.set_title("Chromosome {}".format(chrom), loc='left', weight='bold')
    record.plot(ax=ax)
    pbar.update(n)
pbar.update(10000-n*len(chromosome))
pbar.refresh()
pbar.close()
fig.savefig('{}.svg'.format(out_name), format='svg')