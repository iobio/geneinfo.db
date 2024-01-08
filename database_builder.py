import sqlite3
import csv
import pandas as pd
from data_transformer_refseq import process_refseq_gff
from data_transformer_gencode import process_gencode_gff3
from mane_transcript_update import process_mane_gff
from data_transformer_hgnc import process_hgnc_file
import time
import pandas as pd



# Create SQLite database and table
def create_database(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    sql_cmd = '''
                DROP TABLE IF EXISTS genes;
                DROP TABLE IF EXISTS transcripts;
                DROP TABLE IF EXISTS xref_transcript;
                DROP TABLE IF EXISTS mane_transcripts;

                '''
    c.executescript(sql_cmd)

    sql_cmd_2 = '''
                CREATE TABLE IF NOT EXISTS genes
                (chr text, annotation_source text, feature_type text, start int, end int, score text, strand text, phase text, gene_name text, gene_type text, gene_status text, level int, transcripts text, source text, seqid text, species text, build text, gene_key text);

                CREATE INDEX IF NOT EXISTS gene_name ON genes (gene_name);

                CREATE TABLE IF NOT EXISTS transcripts (chr text, annotation_source text, feature_type text, start int, end int, score text, strand text, phase text, transcript_id text, gene_name text, transcript_type text, transcript_status text, level int, features text, source text, seqid text, build text, species text, gene_key text);

                CREATE INDEX IF NOT EXISTS transcript_id on transcripts (transcript_id);

                CREATE INDEX IF NOT EXISTS start on genes (start);
                
                CREATE INDEX IF NOT EXISTS end on genes (end);

                CREATE TABLE IF NOT EXISTS xref_transcript (gencode_id text, refseq_id text, species text, build text, gene_name text);

                CREATE INDEX IF NOT EXISTS idx_transcripts_composite on transcripts(gene_name,source,species,build);

                CREATE INDEX IF NOT EXISTS idx_xref_transcript_gencode on xref_transcript(gencode_id,species,build);

                CREATE INDEX IF NOT EXISTS idx_xref_transcript_refseq on xref_transcript(refseq_id,species,build);

              '''
    
    c.executescript(sql_cmd_2)

    sql_cmd_3 =   '''
                DROP TABLE IF EXISTS gene_symbol;
                DROP TABLE IF EXISTS gene_alias;
                '''
    c.executescript(sql_cmd_3)

    sql_cmd_4 = '''
                CREATE TABLE IF NOT EXISTS gene_symbol
                 (gene_symbol text NOT NULL,
                  hgnc_id text,
                  gencode_id text,
                  refseq_id text);
                CREATE INDEX IF NOT EXISTS idx_gene_symbol ON gene_symbol (gene_symbol);

                CREATE TABLE IF NOT EXISTS gene_alias
                 (gene_symbol text NOT NULL,
                  alias_symbol text NOT NULL,
                  is_current boolean);
                CREATE INDEX IF NOT EXISTS idx_alias_symbol on gene_alias (alias_symbol);
                CREATE INDEX IF NOT EXISTS idx_gene_symbol on gene_alias (gene_symbol);
                '''

    c.executescript(sql_cmd_4)
    conn.commit()

    sql_cmd_5 = '''
                ALTER TABLE genes ADD COLUMN gene_symbol text;
                ALTER TABLE genes ADD COLUMN alias_symbol text;
                '''
    try:
        c.executescript(sql_cmd_5)
    except:
        print("info: genes table already has columns gene_symbol and alias_symbol.")

    sql_cmd_6 = '''
                CREATE INDEX IF NOT EXISTS idx_alias_symbol on genes (alias_symbol);
                CREATE INDEX IF NOT EXISTS idx_gene_symbol  on genes (gene_symbol);
                '''
    c.executescript(sql_cmd_6)

 
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

#
# Insert hgnc gene symbol data into SQLite database
#
def insert_gene_symbols_into_database(db_name, gene_csv, alias_csv):
    conn = sqlite3.connect(db_name)
    c    = conn.cursor()

    # load the data into a Pandas DataFrame
    gene  = pd.read_csv(gene_csv)
    alias = pd.read_csv(alias_csv)

    # Remove the index column from the DataFrames
    gene  = gene.drop(columns=['Unnamed: 0'])
    alias = alias.drop(columns=['Unnamed: 0'])

    # write the data to a sqlite table
    gene.to_sql( 'gene_symbol',  conn, if_exists='append', index=False)
    alias.to_sql('gene_alias', conn, if_exists='append', index=False)


    conn.commit()
    conn.close()

# Un-escape alias that would have been wrongly interpreted as NULL
# in the pandas.to_sql function
def update_escaped_alias(db_name):
    conn = sqlite3.connect(db_name)
    c    = conn.cursor()

    sql_cmd = "UPDATE gene_alias set alias_symbol = 'NaN' WHERE alias_symbol = '!NaN'"
    c.executescript(sql_cmd)

    sql_cmd = "UPDATE gene_alias set alias_symbol = 'NA' WHERE alias_symbol = '!NA'"
    c.executescript(sql_cmd)

    conn.commit()
    conn.close()

# Update the genes table with cross-reference to gene_symbol and alias_symbol
def update_genes_with_gene_symbol(db_name):
    conn = sqlite3.connect(db_name)
    c    = conn.cursor()

    sql_cmd = "UPDATE genes " + \
              "SET gene_symbol = " + \
              "(SELECT gene_symbol FROM gene_symbol s WHERE genes.gene_name=s.gene_symbol)"

    c.executescript(sql_cmd)

    sql_cmd = "UPDATE genes " + \
              "SET alias_symbol = " + \
              "(SELECT alias_symbol FROM gene_alias a WHERE genes.gene_name=a.alias_symbol)"
    c.executescript(sql_cmd)

    sql_cmd = "UPDATE genes " + \
          " SET gene_symbol = " + \
          "       (select gene_symbol from gene_alias a " + \
          "        WHERE genes.gene_name = a.alias_symbol) " + \
          "WHERE gene_symbol is NULL"
    c.executescript(sql_cmd)

    conn.commit()
    conn.close()

def populate_genes_and_transcripts(db_name):
    file_path_1 = "./data/gencode.v19.annotation.gff3"
    file_path_2 = "./data/gencode.v44.basic.annotation.gff3"

    file_path_3 = './data/GRCh37_latest_genomic.gff'
    file_path_4 = './data/GRCh38_latest_genomic.gff'

    mane_gff_file_path = "data/MANE.GRCh38.v1.0.ensembl_genomic.gff"
    shell_script_path = "extract_mane_transcripts.sh"

    mane_transcrpts = 'data/mane_transcripts.csv'

    print("process gff")
    start = time.time();
    genes_v19_gencode, transcripts_v19_gencode = process_gencode_gff3(file_path_1)
    genes_v44_gencode, transcripts_v44_gencode = process_gencode_gff3(file_path_2)
    genes_GRCh37_refseq, transcripts_GRCh37_refseq = process_refseq_gff(file_path_3)
    genes_GRCh38_refseq, transcripts_GRCh38_refseq = process_refseq_gff(file_path_4)
    # Comment the ^ above lines and uncomment v below lines if you need to re-run without
    # re-creating the intermediate csv files.
    #genes_v19_gencode       = './data/gencode.v19.annotation_gene.csv'
    #transcripts_v19_gencode = './data/gencode.v19.annotation_transcript.csv'
    #genes_v44_gencode       = './data/gencode.v44.basic.annotation_gene.csv'
    #transcripts_v44_gencode = './data/gencode.v44.basic.annotation_transcript.csv'
    #genes_GRCh37_refseq         = './data/GRCh37_latest_genomic_gene.csv'
    #transcripts_GRCh37_refseq   = './data/GRCh37_latest_genomic_transcript.csv'
    #genes_GRCh38_refseq         = './data/GRCh38_latest_genomic_gene.csv'
    #transcripts_GRCh38_refseq   = './data/GRCh38_latest_genomic_transcript.csv'

    end1 = time.time()
    print("\n\nelapsed seconds", round(end1-start,2))
    print("process mane_select")

    process_mane_gff(mane_gff_file_path, shell_script_path)

    end2 = time.time()
    print("\n\nelapsed seconds", round(end2-end1,2))
    print("populate db")

    insert_data_into_database(db_name, genes_v19_gencode, transcripts_v19_gencode)
    insert_data_into_database(db_name, genes_v44_gencode, transcripts_v44_gencode)
    insert_data_into_database(db_name, genes_GRCh37_refseq, transcripts_GRCh37_refseq)
    insert_data_into_database(db_name, genes_GRCh38_refseq, transcripts_GRCh38_refseq)

    end3 = time.time()
    print("\n\nelapsed seconds", round(end3-end2,2))

    update_mane_select_transcripts(db_name, mane_transcrpts)
    end4 = time.time()
    print("\n\nelapsed seconds", round(end4-end3,2))

def populate_gene_symbols(db_name):
    file_path = "./data/hgnc_complete_set.txt"

    print("process hgnc file", flush=True)
    start = time.time()
    gene_symbol, gene_alias = process_hgnc_file(file_path)

    end1 = time.time()
    print("\n\nelapsed seconds", round(end1-start,2), flush=True)

    print("insert gene symbol data into db", flush=True)
    insert_gene_symbols_into_database(db_name, gene_symbol, gene_alias)
    print("update escaped alias", flush=True)
    update_escaped_alias(db_name)
    print("update genes with gene_alias and alias_symbol")
    update_genes_with_gene_symbol(db_name)

    end2 = time.time()
    print("\n\nelapsed seconds", round(end2-end1,2), flush=True)


def main():

    start = time.time()
    print("Building database..")

    db_name = 'gene.iobio.db'

    create_database(db_name)


    #
    # Populate genes and transcripts from
    # gencode and refseq gff3 files
    #
    populate_genes_and_transcripts(db_name)


    #
    # Now populate the gene symbols and aliases
    # from the hgnc file
    #
    populate_gene_symbols(db_name)

    
    print("Database built successfully!")

if __name__ == '__main__':
    main()
