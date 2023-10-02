import pandas as pd
import time


start = time.time()

def print_dups(df):
	grouped = df.groupby('gene_name').filter(lambda x: len(x) > 1).groupby('gene_name')
	print("dup genes count", grouped.sum('count').shape[0], "\n") 
	for name, group in grouped:
		chrs = ""
		group['chr'].fillna('NULL', inplace=True)
		for row_index, row in group.iterrows():
			chrs += row['chr'] + ' '
		print(name, group.shape[0], chrs)

print("\n\n\ngencode GRCh37")
gencode_df = pd.read_csv("data/gencode.v19.annotation_gene.csv")
print_dups(gencode_df)


print("\n\n\ngencode GRCh38")
gencode_df = pd.read_csv("data/gencode.v41.annotation_gene.csv")
print_dups(gencode_df)



print("\n\n\nrefseq GRCh37")
refseq_df = pd.read_csv("data/GRCh37_latest_genomic_gene.csv")
print_dups(refseq_df)


print("\n\n\nrefseq GRCh38")
refseq_df = pd.read_csv("data/GRCh38_latest_genomic_gene.csv")
print_dups(refseq_df)

end = time.time()
print("\n\nelapsed seconds", round(end-start,2))
