#!/usr/bin/env python
# coding: utf-8
import sqlite3
import pandas as pd
import json


# Open a db connection
con = sqlite3.connect("gene.iobio.db")


#
# 
# Create new genes.json based on gene.iobio.db genes
#
#

chromosomes = "('chr1', 'chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY', 'chrM')"

#
# Run a query
#
def run_query(stmt):
    cur = con.cursor()
    res = cur.execute(stmt)
    results = res.fetchall()
    return results

#
# Keep track off all of the warnings encounted for each gene
#
def load_error_gene(error_gene_map, gene_name, source, error):

    error_obj = error_gene_map.get(gene_name, {"gencode": [], "refseq": []})
    error_obj[source].append(error)




#
# Create a gene map based on the sql query that gets
# all distinct entries of gene name, build, and source
#
print("running db query", flush=True)
stmt = "SELECT distinct gene_name, source, build, group_concat(distinct chr)  from genes GROUP BY gene_name, build, source"
results = run_query(stmt)
print("number of rows", len(results), flush=True)


gene_map = {}
gene_names = [];
for row in results:
    gene_name = row[0]
    source = row[1]
    build = row[2]
    chroms = row[3]

    if gene_name in gene_map:
        gene = gene_map[gene_name] 
    else:
        gene = { 
                'GRCh37': {'gencode': 0, 'refseq': 0}, 
                'GRCh38': {'gencode': 0, 'refseq': 0}
                }
        gene_map[gene_name] = gene

    gene[build][source] = len(chroms.split(","))
    
    if gene_name not in gene_names:
        gene_names.append(gene_name)



# Create the genes.json by looping through all of the genes in the
# map.
print("creating genes.json", flush=True)
gene_list = []
gene_names.sort()
for gene_name in gene_names:
    gene_obj = gene_map.get(gene_name, {})


    # if the gene exists in both builds and sources and has no errors,
    # create a simple object with gene_name
    if gene_obj['GRCh37']['gencode'] == 1 and gene_obj['GRCh37']['refseq'] == 1 and gene_obj['GRCh38']['gencode'] == 1 and gene_obj['GRCh38']['refseq'] == 1:
        the_gene_obj = {"gn": gene_name}
        #the_gene_obj = [gene_name]
        #the_gene_obj = gene_name
    # otherwise, create a more complex object that indicates which
    # sources and build have this gene and the errors for each
    # source and build of the gene.
    else:
        the_gene_obj = {"gn": gene_name, "d": []}

        build_info = []
        for build in ['GRCh37', 'GRCh38']:
            source_info = []
            for source in ['gencode', 'refseq']:                
                source_info.append(gene_obj[build][source])
            build_info.append(source_info)
        the_gene_obj["d"] = build_info


    gene_list.append(the_gene_obj)

#
# Write out genes.json
#
print("writing json file", len(gene_list), flush=True)
with open('genes.json', 'w') as fp:
    json.dump(gene_list, fp)    

