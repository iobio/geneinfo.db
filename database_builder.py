import sqlite3
import csv
import pandas as pd
from data_transformer_refseq import process_refseq_gff
from data_transformer_gencode import process_gencode_gff3
from mane_transcript_update import process_mane_gff



# Create SQLite database and table
def create_database(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    sql_cmd = '''
                CREATE TABLE IF NOT EXISTS genes
                (chr text, annotation_source text, feature_type text, start int, end int, score text, strand text, phase text, gene_name text, gene_type text, gene_status text, level int, transcripts text, source text, seqid text, species text, build text);

                CREATE INDEX IF NOT EXISTS gene_name ON genes (gene_name);

                CREATE TABLE IF NOT EXISTS transcripts (chr text, annotation_source text, feature_type text, start int, end int, score text, strand text, phase text, transcript_id text, gene_name text, transcript_type text, transcript_status text, level int, features text, source text, seqid text, build text, species text);

                CREATE INDEX IF NOT EXISTS transcript_id on transcripts (transcript_id);

                CREATE INDEX IF NOT EXISTS start on genes (start);
                
                CREATE INDEX IF NOT EXISTS end on genes (end);

                CREATE TABLE IF NOT EXISTS xref_transcript (gencode_id text, refseq_id text, species text, build text, gene_name text);

                CREATE INDEX IF NOT EXISTS idx_transcripts_composite on transcripts(gene_name,source,species,build);

                CREATE INDEX IF NOT EXISTS idx_xref_transcript_gencode on xref_transcript(gencode_id,species,build);

                CREATE INDEX IF NOT EXISTS idx_xref_transcript_refseq on xref_transcript(refseq_id,species,build);

              '''
    
    c.executescript(sql_cmd)
 
    conn.commit()
    conn.close()


# Update the is_mane_select column in the transcripts table
def update_mane_select_transcripts(db_name, mane_transcripts_csv):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    sql_cmd = '''
                CREATE TABLE IF NOT EXISTS mane_transcripts (gene_name text, transcript_id text, transcript_id_refseq text);
                CREATE INDEX idx_mane_transcripts on mane_transcripts(transcript_id);
                CREATE INDEX idx_mane_transcripts_refseq on mane_transcripts(transcript_id_refseq);

              '''
    
    c.executescript(sql_cmd)

    mane_transcripts_df = pd.read_csv(mane_transcripts_csv)  
    mane_transcripts_df.to_sql('mane_transcripts', conn, if_exists='append', index=False)
    
    sql_cmd_2 = '''
                ALTER TABLE transcripts ADD COLUMN is_mane_select TEXT;

                UPDATE transcripts
                SET is_mane_select = 'true'
                WHERE EXISTS (
                    SELECT 1
                    FROM mane_transcripts
                    WHERE (mane_transcripts.transcript_id = transcripts.transcript_id AND transcripts.build = 'GRCh38' AND transcripts.source = 'gencode')
                    OR (mane_transcripts.transcript_id_refseq = transcripts.transcript_id AND transcripts.build = 'GRCh38' AND transcripts.source = 'refseq')
                );

                DROP TABLE mane_transcripts;
                
                '''
 
    c.executescript(sql_cmd_2)

    conn.commit()
    conn.close()


# Insert data into SQLite database
def insert_data_into_database(db_name, genes_csv, transcripts_csv):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    # load the data into a Pandas DataFrame
    genes = pd.read_csv(genes_csv)
    transcripts = pd.read_csv(transcripts_csv)

    # Remove the index column from the DataFrames
    genes = genes.drop(columns=['Unnamed: 0'])
    transcripts = transcripts.drop(columns=['Unnamed: 0'])
 
    # write the data to a sqlite table
    genes.to_sql('genes', conn, if_exists='append', index = False)
    transcripts.to_sql('transcripts', conn, if_exists='append', index=False)
    conn.commit()
    conn.close()


def main():

    print("Building database..")

    # genes_v19_gencode = 'data/gencode.v19.annotation_gene.csv'
    # transcripts_v19_gencode = 'data/gencode.v19.annotation_transcript.csv'

    # genes_v41_gencode = 'data/gencode.v41.annotation_gene.csv'
    # transcripts_v41_gencode = 'data/gencode.v41.annotation_transcript.csv'

    # genes_GRCh37_refseq = 'data/GRCh37_latest_genomic_gene.csv'
    # transcripts_GRCh37_refseq = 'data/GRCh37_latest_genomic_transcript.csv'

    # genes_GRCh38_refseq = 'data/GRCh38_latest_genomic_gene.csv'
    # transcripts_GRCh38_refseq = 'data/GRCh38_latest_genomic_transcript.csv'

    mane_transcrpts = 'data/mane_transcripts.csv'

    db_name = 'gene.iobio.db'

    file_path_1 = "./data/gencode.v19.annotation.gff3"
    file_path_2 = "./data/gencode.v41.annotation.gff3"

    file_path_3 = './data/GRCh37_latest_genomic.gff'
    file_path_4 = './data/GRCh38_latest_genomic.gff'

    mane_gff_file_path = "data/MANE.GRCh38.v1.0.ensembl_genomic.gff"
    shell_script_path = "extract_mane_transcripts.sh"

    genes_v19_gencode, transcripts_v19_gencode = process_gencode_gff3(file_path_1)
    genes_v41_gencode, transcripts_v41_gencode = process_gencode_gff3(file_path_2)
    genes_GRCh37_refseq, transcripts_GRCh37_refseq = process_refseq_gff(file_path_3)
    genes_GRCh38_refseq, transcripts_GRCh38_refseq = process_refseq_gff(file_path_4)
    process_mane_gff(mane_gff_file_path, shell_script_path)
  
    create_database(db_name)
    insert_data_into_database(db_name, genes_v19_gencode, transcripts_v19_gencode)
    insert_data_into_database(db_name, genes_v41_gencode, transcripts_v41_gencode)
    insert_data_into_database(db_name, genes_GRCh37_refseq, transcripts_GRCh37_refseq)
    insert_data_into_database(db_name, genes_GRCh38_refseq, transcripts_GRCh38_refseq)
    update_mane_select_transcripts(db_name, mane_transcrpts)
    
    print("Database built successfully!")

if __name__ == '__main__':
    main()
