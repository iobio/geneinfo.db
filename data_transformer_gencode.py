#!/usr/bin/env python
# coding: utf-8

import pandas as pd 
import gffpandas.gffpandas as gffpd
import os

def process_gencode_gff3(file_path):

    base_file_name = os.path.splitext(file_path)[0]
    output_file_gene = f"{base_file_name}_gene.csv"
    output_file_transcript = f"{base_file_name}_transcript.csv"

    annotation = gffpd.read_gff3(file_path)
    attr_to_columns = annotation.attributes_to_columns()

    # Filtered out the rows with gene type
    filtered_gene_df = attr_to_columns[attr_to_columns["type"] == 'gene']

    # Filtered out the rows with transcript type
    filtered_transcript_df = attr_to_columns[attr_to_columns["type"] == 'transcript']

    # Create a new column 'transcripts' 
    filtered_gene_df.loc[:,'transcripts'] = None

    # Mapping gene_name to a list of transcript_ids
    transcript_dict = filtered_transcript_df.groupby('gene_name')['transcript_id'].apply(list).to_dict()

    filtered_gene_df.loc[:,'transcripts'] = filtered_gene_df['gene_name'].map(transcript_dict)

    filtered_gene_df.loc[:,'chr'] = filtered_gene_df['seq_id']

    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid'}
    filtered_gene_df.rename(columns=new_column_names, inplace=True)

    # Created other new Columns
    filtered_gene_df.loc[:,'source'] = "gencode"
    filtered_gene_df.loc[:,'species'] = "homo_sapiens"

    if "GRCh38" in annotation.header:
        filtered_gene_df.loc[:,'gene_status'] = "."
        
    if "GRCh37" in annotation.header:
        filtered_gene_df.loc[:,'build'] = "GRCh37"
    else:
        filtered_gene_df.loc[:,'build'] = "GRCh38"
        

    # List of columns to keep
    columns_to_keep = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 
                    'gene_name', 'gene_type', 'gene_status', 'level', 'transcripts', 'source', 'seqid', 
                    'species', 'build']

    filtered_gene_df_final = filtered_gene_df[columns_to_keep]
    filtered_gene_df_final.to_csv(output_file_gene)
    print("Writing gene file done!")

    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid'}
    filtered_transcript_df.rename(columns=new_column_names, inplace=True)

    filtered_transcript_df.loc[:,'chr'] = filtered_transcript_df['seqid']

    # Created other new Columns
    filtered_transcript_df.loc[:,'source'] = "gencode"
    filtered_transcript_df.loc[:,'species'] = "homo_sapiens"

    # Create a new column 'features' 
    filtered_transcript_df.loc[:,'features'] = None

    if "GRCh38" in annotation.header:
        filtered_transcript_df['transcript_status'] = "."
        
    if "GRCh37" in annotation.header:
        filtered_transcript_df.loc[:,'build'] = "GRCh37"
    else:
        filtered_transcript_df.loc[:,'build'] = "GRCh38"

    # List of columns to keep
    columns_to_keep = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 
                    'transcript_id', 'gene_name', 'transcript_type', 'transcript_status', 'level', 'features',
                    'source', 'seqid', 'build', 'species']

    filtered_transcript_df = filtered_transcript_df[columns_to_keep]

    #Filter out the rows by type column with values: 'exon','CDS','stop_codon','start_codon' 'UTR' 
    feature_type_to_keep = [ 'exon','CDS','stop_codon','start_codon','UTR']

    filtered_feature_type_df = attr_to_columns[attr_to_columns['type'].isin(feature_type_to_keep)] 

    # filtered_feature_type_df['chr'] = filtered_feature_type_df['seq_id'].copy()
    filtered_feature_type_df.loc[:,'chr'] = filtered_feature_type_df['seq_id']



    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid'}
    filtered_feature_type_df.rename(columns=new_column_names, inplace=True)

    # List of columns to keep
    columns_to_keep = ["chr","seqid","annotation_source","feature_type","start","end","score",
                    "strand","phase","transcript_id"]

    filtered_feature_type_df = filtered_feature_type_df[columns_to_keep]

    feature_columns = ["chr", "seqid", "annotation_source", "feature_type", "start", "end", "score",
                    "strand", "phase", "transcript_id"]

    # Group by 'transcript_id' and create a dictionary of features for each transcript
    features_dict = filtered_feature_type_df.groupby('transcript_id')[feature_columns].apply(lambda x: x.to_dict('records')).to_dict()

    # Map the 'features' column using the 'transcript_id'
    filtered_transcript_df['features'] = filtered_transcript_df['transcript_id'].map(features_dict)

    filtered_transcript_df.to_csv(output_file_transcript)

    print("Writing transcript file done!")

    return output_file_gene, output_file_transcript





