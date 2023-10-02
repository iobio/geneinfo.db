#!/usr/bin/env python
# coding: utf-8

import json
import pandas as pd 
import gffpandas.gffpandas as gffpd
import os
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
    print("refseq ", table_name, build)
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





def process_refseq_gff(file_path):

    base_file_name = os.path.splitext(file_path)[0]
    output_file_gene = f"{base_file_name}_gene.csv"
    output_file_transcript = f"{base_file_name}_transcript.csv"

    print("reading gff")
    annotation = gffpd.read_gff3(file_path)
    attr_to_columns = annotation.attributes_to_columns()

    # Filtered out the rows with gene type
    filtered_gene_df = attr_to_columns[attr_to_columns['type'] == 'gene'].copy()

    # Under type column, keep certain types
    selected_types = ['transcript','mRNA','miRNA','ncRNA','rRNA','snoRNA','snRNA','tRNA','misc_RNA',
                        'processed_transcript','primary-transcript','scRNA']

    # Create a new DataFrame based on the selected values
    filtered_transcript_df = attr_to_columns[attr_to_columns['type'].isin(selected_types)].copy()

    # Rename gene_id column to gene_key on the genes dataframe and the transcripts 
    # dataframe. We will use this gene_key to guarantee uniqueness since we can 
    # have mulitiple gene records in the gff with the same gene name
    filtered_gene_df.rename(columns={'ID': 'gene_key'}, inplace=True)
    filtered_transcript_df.rename(columns={'Parent': 'gene_key'}, inplace=True)


    # Create a new column 'transcripts' 
    filtered_gene_df['transcripts'] = None

    # Mapping gene_name to a list of transcript_ids
    transcript_dict = filtered_transcript_df.groupby('gene_key')['transcript_id'].unique().apply(list).to_dict()

    # # Create a new column 'transcripts' by applying a function that maps values from transcript_dict
    # filtered_gene_df['transcripts'] = filtered_gene_df.apply(lambda row: transcript_dict.get((row['gene'], row['seq_id']), []), axis=1)


    filtered_gene_df['transcripts'] = filtered_gene_df['gene_key'].map(transcript_dict)

    # Convert the list of transcripts to a json string
    filtered_gene_df['transcripts'] = filtered_gene_df['transcripts'].apply(lambda x: json.dumps(x))

    filtered_gene_df['chr'] = filtered_gene_df['seq_id']

    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid', 'gene':'gene_name'}
    filtered_gene_df.rename(columns=new_column_names, inplace=True)

    # Created other new Columns
    filtered_gene_df['source'] = "refseq"
    filtered_gene_df['species'] = "homo_sapiens"
    filtered_gene_df['gene_type'] = "."
    filtered_gene_df['gene_status'] = "."
    filtered_gene_df['level'] = "."
        
    if "GRCh37" in annotation.header:
        filtered_gene_df['build'] = "GRCh37"
    else:
        filtered_gene_df['build'] = "GRCh38"
        

    # List of columns to keep
    columns_to_keep = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 
                    'gene_name', 'gene_type', 'gene_status', 'level', 'transcripts', 'source', 'seqid', 
                    'species', 'build', 'gene_key']

    filtered_gene_df_final = filtered_gene_df[columns_to_keep].copy()


    # This is for mapping sequence IDs to chromosomes
    if "GRCh37" in base_file_name:
        file_path_ = 'data/GCF_000001405.25_GRCh37.p13_assembly_report.txt'
    else:
        file_path_ = 'data/GCF_000001405.40_GRCh38.p14_assembly_report.txt'
    columns_to_read = [1, 6, 9]  # Index of the 1st, 7th and 10th columns 
    mapping_df = pd.read_csv(file_path_, comment='#', delimiter='\t', usecols=columns_to_read, header=None)

    specific_values_to_exclude = ['fix-patch', 'novel-patch']

    mapping_df = mapping_df[~mapping_df[1].isin(specific_values_to_exclude)]

    mapping_df = mapping_df.loc[:, [6,9]]

    column_names = ['seqid', 'chr']
    mapping_df.columns = column_names

    seqid_to_chr = dict(zip(mapping_df['seqid'], mapping_df['chr']))

    filtered_gene_df_final['chr'] = filtered_gene_df_final['seqid'].map(seqid_to_chr)

    # Filter out duplicate genes (on non-normal chromosomes)
    filtered_gene_df_final = filter_dup_genes(filtered_gene_df_final, 
                                              'genes', 
                                              'GRCh37' if 'GRCh37' in annotation.header else 'GRCh38')

    start = time.time()
    filtered_gene_df_final.to_csv(output_file_gene)
    end = time.time()
    print(str(base_file_name) + "_gene file writing done!")
    print("time to write csv file", round(end-start,2), "seconds")


    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid','gene':'gene_name'}
    filtered_transcript_df.rename(columns=new_column_names, inplace=True)

    # filtered_transcript_df['chr'] = filtered_transcript_df['seqid']
    filtered_transcript_df['chr'] = filtered_transcript_df['seqid'].map(seqid_to_chr)

    # Filter out duplicate genes (on non-normal chromosomes)
    filtered_transcript_df = filter_dup_genes(filtered_transcript_df,
                                              'transcripts', 
                                              'GRCh37' if 'GRCh37' in annotation.header else 'GRCh38')

    # Created other new Columns
    filtered_transcript_df['source'] = "refseq"
    filtered_transcript_df['species'] = "homo_sapiens"
    filtered_transcript_df['transcript_type'] = "."
    filtered_transcript_df['transcript_status'] = "."
    filtered_transcript_df['level'] = "."

    # Create a new column 'features' 
    filtered_transcript_df['features'] = None
        
    if "GRCh37" in annotation.header:
        filtered_transcript_df['build'] = "GRCh37"
    else:
        filtered_transcript_df['build'] = "GRCh38"

    # List of columns to keep
    columns_to_keep = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 
                    'transcript_id', 'gene_name', 'transcript_type', 'transcript_status', 'level', 'features',
                    'source', 'seqid', 'build', 'species', 'gene_key']


    filtered_transcript_df = filtered_transcript_df[columns_to_keep]


    #Filter out the rows by type column with values: 'exon','CDS','stop_codon','start_codon' 'UTR' 
    feature_type_to_keep = [ 'exon','CDS','stop_codon','start_codon','UTR']

    filtered_feature_type_df = attr_to_columns[attr_to_columns['type'].isin(feature_type_to_keep)].copy()

    filtered_feature_type_df['chr'] = filtered_feature_type_df['seq_id']

    filtered_feature_type_df['transcript_id'] = filtered_feature_type_df['Parent'].str.replace('gene-', '').str.replace('rna-', '')

    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid'}
    filtered_feature_type_df.rename(columns=new_column_names, inplace=True)

    # List of columns to keep
    columns_to_keep = ["chr","seqid","annotation_source","feature_type","start","end","score",
                    "strand","phase","transcript_id"]

    filtered_feature_type_df = filtered_feature_type_df[columns_to_keep]
    filtered_feature_type_df['chr'] = filtered_feature_type_df['seqid'].map(seqid_to_chr)


    # Group data by transcript_id
    grouped = filtered_feature_type_df.groupby("transcript_id")

    # Calculate UTR for each transcript
    utr_data = []
    for transcript_id, group_data in grouped:
        exon_data = group_data[group_data["feature_type"] == "exon"]
        cds_data = group_data[group_data["feature_type"] == "CDS"]
        
        # Calculate UTR regions
        for _, exon_row in exon_data.iterrows():
            exon_start = exon_row["start"]
            exon_end = exon_row["end"]
            has_cds = False
            
            # Check if exon has corresponding CDS regions
            for _, cds_row in cds_data.iterrows():
                cds_start = cds_row["start"]
                cds_end = cds_row["end"]
                if exon_start <= cds_start and exon_end >= cds_end:
                    has_cds = True
                    utr_5_start = exon_start
                    utr_5_end = cds_start
                    utr_3_start = cds_end
                    utr_3_end = exon_end
                    if utr_5_start != utr_5_end:
                        utr_data.append((exon_row["chr"],exon_row["seqid"],exon_row["annotation_source"],"UTR", utr_5_start, utr_5_end,exon_row["strand"],exon_row["score"],".",transcript_id))
                    if utr_3_start != utr_3_end:
                        utr_data.append((exon_row["chr"],exon_row["seqid"],exon_row["annotation_source"],"UTR", utr_3_start, utr_3_end,exon_row["strand"],exon_row["score"],".",transcript_id))
                    break
                    
            if not has_cds:
                if exon_start != exon_end:
                    utr_data.append((exon_row["chr"],exon_row["seqid"],exon_row["annotation_source"],"UTR", exon_start, exon_end,exon_row["strand"],exon_row["score"],".",transcript_id))


    utr_df = pd.DataFrame(utr_data, columns=columns_to_keep)

    # Concatenate UTR, exon and CDS
    filtered_feature_utr_df = pd.concat([filtered_feature_type_df, utr_df],ignore_index=True)


    feature_columns = ["chr", "seqid", "annotation_source", "feature_type", "start", "end", "score",
                    "strand", "phase", "transcript_id"]

    # Group by 'transcript_id' and create a dictionary of features for each transcript
    features_dict = filtered_feature_utr_df.groupby('transcript_id')[feature_columns].apply(lambda x: x.to_dict('records')).to_dict()
    
    # Map the 'features' column using the 'transcript_id'
    filtered_transcript_df['features'] = filtered_transcript_df['transcript_id'].map(features_dict)

    # Convert the list of features to a json string
    filtered_transcript_df['features'] = filtered_transcript_df['features'].apply(lambda x: json.dumps(x))
   
    filtered_transcript_df['feature_type'] = "transcript"

    start = time.time()
    filtered_transcript_df.to_csv(output_file_transcript)
    end = time.time()

    print(str(base_file_name) + "_transcript file writing done!")
    print("time to write csv file", round(end-start,2), "seconds")

    return output_file_gene, output_file_transcript
