#!/usr/bin/env python
# coding: utf-8

import pandas as pd 
import os
import time
import sqlite3


#
# Run a query
#
def run_query(con, stmt):
    cur = con.cursor()
    res = cur.execute(stmt)
    results = res.fetchall()
    return results

# None and NA are treated as NULL in pandas.to_sql. These are actual
# alias names, so we will prepend with a exclamation mark and update
# post populate
def escapeAlias(alias):
    if alias == 'NaN' or alias == 'NA':
        return "!" + alias;
    else:
        return alias;

def process_hgnc_file(file_path):

    # Open a db connection to gene.iobio.db
    con = sqlite3.connect("gene.iobio.db")

    base_file_name    = os.path.splitext(file_path)[0]
    output_file_gene  = f"{base_file_name}_gene_symbol.csv"
    output_file_alias = f"{base_file_name}_gene_alias.csv"

    print('reading hgnc file')
    hgnc = pd.read_table(file_path, sep="\t")

    gene_records = []
    alias_records = []

    for idx, rec in hgnc.iterrows():
        gene_records.append({
            'gene_symbol': rec['symbol'],
            'hgnc_id': rec['hgnc_id'],
            'gencode_id': rec['ensembl_gene_id'],
            'refseq_id': rec['refseq_accession']
        })

        # parse the alias field
        aliases = []
        if type(rec['alias_symbol']) is str:
            aliases = rec['alias_symbol'].split("|")
            for alias in aliases:
                alias_records.append({
                    'gene_symbol': rec['symbol'],
                    'alias_symbol': escapeAlias(alias),
                    'is_current': True
                    })
        prev_aliases = []
        if type(rec['prev_symbol']) is str:
            prev_aliases = rec['prev_symbol'].split("|")
            for alias in prev_aliases:
                alias_records.append({
                    'gene_symbol': rec['symbol'],
                    'alias_symbol': escapeAlias(alias),
                    'is_current': False
                    })

    start = time.time()
    gene_df = pd.DataFrame.from_records(gene_records)
    gene_df.to_csv(output_file_gene)
    end = time.time()
    print(str(base_file_name) + "_gene_symbol file writing done.")

    alias_df = pd.DataFrame.from_records(alias_records)
    alias_df.to_csv(output_file_alias)
    end1 = time.time()
    print(str(base_file_name) + "_gene_alias file writing done.")
    print("time to write csv file", round(end1-end,2), "seconds")


    return output_file_gene, output_file_alias





