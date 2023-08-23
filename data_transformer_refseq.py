#!/usr/bin/env python
# coding: utf-8

import pandas as pd 
import gffpandas.gffpandas as gffpd
import os

def process_refseq_gff(file_path):

    base_file_name = os.path.splitext(file_path)[0]
    output_file_gene = f"{base_file_name}_gene.csv"
    output_file_transcript = f"{base_file_name}_transcript.csv"

    annotation = gffpd.read_gff3(file_path)
    attr_to_columns = annotation.attributes_to_columns()

    #Filtered out the rows with gene type
    filtered_gene_df = attr_to_columns[attr_to_columns['type'] == 'gene'].copy()


    # Under type column, keep certain types
    selected_types = ['transcript','mRNA','miRNA','ncRNA','rRNA','snoRNA','snRNA','tRNA','misc_RNA',
                        'processed_transcript','primary-transcript','scRNA']

    # Create a new DataFrame based on the selected values
    filtered_transcript_df = attr_to_columns[attr_to_columns['type'].isin(selected_types)].copy()

    # Create a new column 'transcripts' 
    filtered_gene_df['transcripts'] = None

    # Mapping gene_name to a list of transcript_ids
    transcript_dict = filtered_transcript_df.groupby('gene')['transcript_id'].unique().apply(list).to_dict()

    # # Create a new column 'transcripts' by applying a function that maps values from transcript_dict
    # filtered_gene_df['transcripts'] = filtered_gene_df.apply(lambda row: transcript_dict.get((row['gene'], row['seq_id']), []), axis=1)


    filtered_gene_df['transcripts'] = filtered_gene_df['gene'].map(transcript_dict)

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
        filtered_gene_df['build'] = "GRCH37"
    else:
        filtered_gene_df['build'] = "GRCH38"
        

    # List of columns to keep
    columns_to_keep = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 
                    'gene_name', 'gene_type', 'gene_status', 'level', 'transcripts', 'source', 'seqid', 
                    'species', 'build']

    filtered_gene_df_final = filtered_gene_df[columns_to_keep].copy()


    # This is for mapping sequence IDs to chromosomes
    if "GRCh37" in base_file_name:
        file_path_ = 'data/GCF_000001405.25_GRCh37.p13_assembly_report.txt'
    else:
        file_path_ = 'data/GCF_000001405.40_GRCh38.p14_assembly_report.txt'
    columns_to_read = [6, 9]  # Index of the 7th and 10th columns 
    mapping_df = pd.read_csv(file_path_, comment='#', delimiter='\t', usecols=columns_to_read, header=None)

    column_names = ['seqid', 'chr']
    mapping_df.columns = column_names

    seqid_to_chr = dict(zip(mapping_df['seqid'], mapping_df['chr']))

    filtered_gene_df_final['chr'] = filtered_gene_df_final['seqid'].map(seqid_to_chr)

    filtered_gene_df_final.to_csv(output_file_gene)
    print("Writing gene file done!")


    # Rename columns
    new_column_names = {'source': 'annotation_source', 'type': 'feature_type', 'seq_id':'seqid','gene':'gene_name'}
    filtered_transcript_df.rename(columns=new_column_names, inplace=True)

    # filtered_transcript_df['chr'] = filtered_transcript_df['seqid']
    filtered_transcript_df['chr'] = filtered_transcript_df['seqid'].map(seqid_to_chr)


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
                    'source', 'seqid', 'build', 'species']


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

    filtered_transcript_df['feature_type'] = "transcript"

    filtered_transcript_df.to_csv(output_file_transcript)

    print("Writing transcript file done!")

    return output_file_gene, output_file_transcript