#!/usr/bin/env python
# coding: utf-8
import sqlite3
import pandas as pd
import json

#
#
# Create new genes.json based on gene.lookup.db genes
#
#

#
# Run a query
#
def run_query(conn, stmt):
    cur = conn.cursor()
    res = cur.execute(stmt)
    results = res.fetchall()
    return results



# Open a db connection
conn1 = sqlite3.connect("gene.iobio.db")

#
# Create a gene map based on the sql query that gets
# all distinct entries of gene symbol, build, and source
# as well as all distinct entries for alias_symbol,
# build, and source
#
stmt =  '''
        SELECT
            g.gene_name,
            g.build,
            g.source,
            LENGTH(g.transcripts) as transcript_count,
            gs.gene_symbol,
            GROUP_CONCAT(ga.alias_symbol) AS aliases
        FROM genes g
        LEFT JOIN gene_symbol gs
          ON g.gene_symbol = gs.gene_symbol
        LEFT JOIN gene_alias ga
          ON gs.gene_symbol = ga.gene_symbol and ga.alias_symbol != g.gene_name
        GROUP BY g.gene_name, g.source, g.build;
         '''

results = run_query(conn1, stmt)


gene_map = {}
gene_names = [];
count = 0
for row in results:
    gene_name        = row[0]
    build            = row[1]
    source           = row[2]
    transcript_count = row[3]
    gene_symbol      = row[4]
    aliases          = row[5]

    if gene_name is None:
        print("Warning, invalid gene name", row)

    if gene_name in gene_map:
        gene = gene_map[gene_name]
    else:
        gene = {
                'GRCh37': {'gencode': 0, 'refseq': 0},
                'GRCh38': {'gencode': 0, 'refseq': 0},
                }
        gene_map[gene_name] = gene

    if build is not None and source is not None and transcript_count is not None:
        gene[build][source] = transcript_count

    # Add the gene_symbol to the aliases (if the gene_symbol is
    # different than the gene name)
    if gene_symbol is not None and gene_symbol != gene_name:
        if aliases is None:
            aliases = gene_symbol
        elif gene_symbol not in aliases:
            aliases = gene_symbol + "," + aliases

    if aliases is not None and len(aliases) > 0:
        gene['aliases'] = aliases

    if gene_name is not None and gene_name not in gene_names:
        gene_names.append(gene_name)

    #if len(gene_names) > 10000:
    #    break

    count += 1





# Create the genes.json by looping through all of the genes in the
# map.
print("creating genes.json", flush=True)
gene_list = []
gene_names.sort()
for gene_name in gene_names:
    gene_obj = gene_map.get(gene_name, {})


    # create an object that indicates which
    # sources and build have this gene
    the_gene_obj = {"gn": gene_name, "d": []}
    if 'aliases' in gene_obj:
        the_gene_obj['a'] = gene_obj['aliases']

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

