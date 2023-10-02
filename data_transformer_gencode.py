#!/usr/bin/env python
# coding: utf-8

import pandas as pd 
import gffpandas.gffpandas as gffpd
import os
import json
import time

normal_chrs = ['chr1',
 'chr2',
 'chr3',
 'chr4',
 'chr5',
 'chr6',
 'chr7',
 'chr8',
 'chr9',
 'chr10',
 'chr11',
 'chr12',
 'chr13',
 'chr14',
 'chr15',
 'chr16',
 'chr17',
 'chr18',
 'chr19',
 'chr20',
 'chr21',
 'chr22',
 'chrX',
 'chrY',
 'chrM'
  ]

def filter_dup_genes(df, table_name, build):
    start = time.time()

    cols = df.columns
    print("************************************")
    print("gencode ", table_name, build)
    print("************************************")
    count0 = df.shape[0]
    print("rows                  ", count0)
    df = df[~df['chr'].isna()]
    count1 = df.shape[0]
    print("after null chr removed", count1, '   ', count0-count1, 'removed')

    to_discard_all= []
    grouped = df.groupby('gene_name')
    grouped_and_filtered = grouped.filter(lambda x: len(x) > 1).groupby('gene_name')
    for name, group in grouped_and_filtered:
        to_discard = []
        to_keep = []
        for row_idx, row in group.iterrows():
            if row['chr'] in normal_chrs:
                to_keep.append(row['gene_key'])
            else:
                to_discard.append(row['gene_key'])
        if len(to_keep) > 0:
            for key in to_discard:
                to_discard_all.append(key)
    print("keys to remove        ", len(to_discard_all))
    df_filtered = df[~df['gene_key'].isin(to_discard_all)]
    count2 = df_filtered.shape[0]
    print("after filter dup genes", count2, '   ', count1-count2, 'removed')

    end = time.time()
    print("\nfilter dup genes", round(end-start,2), "seconds")
    
    return df_filtered[cols]





def process_gencode_gff3(file_path):

    base_file_name = os.path.splitext(file_path)[0]
    output_file_gene = f"{base_file_name}_gene.csv"
    output_file_transcript = f"{base_file_name}_transcript.csv"

    annotation = gffpd.read_gff3(file_path)
    attr_to_columns = annotation.attributes_to_columns()

    # Filtered out the rows with gene type
    filtered_gene_df = attr_to_columns[attr_to_columns["type"] == 'gene'].copy()

    # Filtered out the rows with transcript type
    filtered_transcript_df = attr_to_columns[attr_to_columns["type"] == 'transcript'].copy()

    # Create a new column 'transcripts' 
    filtered_gene_df['transcripts'] = None

    # Rename gene_id column to gene_key on the genes dataframe and the transcripts 
    # dataframe. We will use this gene_key to guarantee uniqueness since we can 
    # have mulitiple gene records in the gff with the same gene name
    new_column_names = {'gene_id': 'gene_key'}
    filtered_gene_df.rename(columns=new_column_names, inplace=True)
    filtered_transcript_df.rename(columns=new_column_names, inplace=True)

    # Mapping gene_name to a list of transcript_ids
    transcript_dict = filtered_transcript_df.groupby('gene_key')['transcript_id'].apply(list).to_dict()

    filtered_gene_df['transcripts'] = filtered_gene_df['gene_key'].map(transcript_dict)

    # Convert the list of transcripts to a json string
    filtered_gene_df['transcripts'] = filtered_gene_df['transcripts'].apply(lambda x: json.dumps(x))
    
    filtered_gene_df['chr'] = filtered_gene_df['seq_id']

    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid'}
    filtered_gene_df.rename(columns=new_column_names, inplace=True)

    # Created other new Columns
    filtered_gene_df['source'] = "gencode"
    filtered_gene_df['species'] = "homo_sapiens"

    if "GRCh38" in annotation.header:
        filtered_gene_df['gene_status'] = "."
        
    if "GRCh37" in annotation.header:
        filtered_gene_df['build'] = "GRCh37"
    else:
        filtered_gene_df['build'] = "GRCh38"
        

    # List of columns to keep
    columns_to_keep = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 
                    'gene_name', 'gene_type', 'gene_status', 'level', 'transcripts', 'source', 'seqid', 
                    'species', 'build', 'gene_key']

    filtered_gene_df_final = filtered_gene_df[columns_to_keep]

    # Filter out dup genes (on non-regular chromosomes)
    filtered_gene_df_final = filter_dup_genes(filtered_gene_df_final, 
                                              'genes',
                                              'GRCh37' if 'GRCh37' in annotation.header else 'GRCh38')

    start = time.time()
    filtered_gene_df_final.to_csv(output_file_gene)
    end = time.time()
    print(str(base_file_name) + "_gene file writing done!")
    print("time to write csv file", round(end-start,2), "seconds")

    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid', 'gene_id': 'gene_key'}
    filtered_transcript_df.rename(columns=new_column_names, inplace=True)

    filtered_transcript_df['chr'] = filtered_transcript_df['seqid']

    # Created other new Columns
    filtered_transcript_df['source'] = "gencode"
    filtered_transcript_df['species'] = "homo_sapiens"

    # Create a new column 'features' 
    filtered_transcript_df['features'] = None

    if "GRCh38" in annotation.header:
        filtered_transcript_df['transcript_status'] = "."
        
    if "GRCh37" in annotation.header:
        filtered_transcript_df['build'] = "GRCh37"
    else:
        filtered_transcript_df['build'] = "GRCh38"

    # List of columns to keep
    columns_to_keep = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 
                    'transcript_id', 'gene_name', 'transcript_type', 'transcript_status', 'level', 'features',
                    'source', 'seqid', 'build', 'species', 'gene_key']

    filtered_transcript_df = filtered_transcript_df[columns_to_keep].copy()

    #Filter out the rows by type column with values: 'exon','CDS','stop_codon','start_codon' 'UTR' 
    feature_type_to_keep = [ 'exon','CDS','stop_codon','start_codon','UTR','three_prime_UTR','five_prime_UTR']

    filtered_feature_type_df = attr_to_columns[attr_to_columns['type'].isin(feature_type_to_keep)].copy()
    filtered_feature_type_df.loc[filtered_feature_type_df['type'] == 'three_prime_UTR', 'type'] = 'UTR'
    filtered_feature_type_df.loc[filtered_feature_type_df['type'] == 'five_prime_UTR', 'type'] = 'UTR'
   

    filtered_feature_type_df['chr'] = filtered_feature_type_df['seq_id']

    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid', 'gene_id': 'gene_key'}
    filtered_feature_type_df.rename(columns=new_column_names, inplace=True)

    # List of columns to keep
    columns_to_keep = ["chr","seqid","annotation_source","feature_type","start","end","score",
                    "strand","phase","transcript_id", "gene_key"]

    filtered_feature_type_df = filtered_feature_type_df[columns_to_keep]

    feature_columns = ["chr", "seqid", "annotation_source", "feature_type", "start", "end", "score",
                    "strand", "phase", "transcript_id"]

    # Group by 'transcript_id' and create a dictionary of features for each transcript
    features_dict = filtered_feature_type_df.groupby('transcript_id')[feature_columns].apply(lambda x: x.to_dict('records')).to_dict()

    # Map the 'features' column using the 'transcript_id'
    filtered_transcript_df['features'] = filtered_transcript_df['transcript_id'].map(features_dict)

    # Convert the list of features to a json string
    filtered_transcript_df['features'] = filtered_transcript_df['features'].apply(lambda x: json.dumps(x))

    # Filter out dup genes (on non-regular chromosomes)
    filtered_transcript_df = filter_dup_genes(filtered_transcript_df,
                                              'transcripts',
                                              'GRCh37' if 'GRCh37' in annotation.header else 'GRCh38')

    start = time.time()
    filtered_transcript_df.to_csv(output_file_transcript)
    end = time.time()
    print(str(base_file_name) + "_transcript file writing done!")
    print("time to write csv file", round(end-start,2), "seconds")


    return output_file_gene, output_file_transcript





