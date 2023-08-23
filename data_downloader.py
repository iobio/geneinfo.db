import ftplib
import os
import gzip

def download_and_extract_files(ftp_server, remote_paths, local_folder):
    ftp = ftplib.FTP(ftp_server)
    ftp.login()

    for remote_path in remote_paths:
        filename = os.path.basename(remote_path)
        local_path = os.path.join(local_folder, filename)

        # Download file
        with open(local_path, 'wb') as local_file:
            ftp.retrbinary('RETR ' + remote_path, local_file.write)
        
        if remote_path.endswith('.gz'):

            # Extract the Gzip file contents
            extracted_path = os.path.join(local_folder, os.path.splitext(filename)[0])
            
            with gzip.open(local_path, 'rb') as gzip_file:
                extracted_content = gzip_file.read()

            # Write the extracted file
            with open(extracted_path, 'wb') as extracted_file:
                extracted_file.write(extracted_content)

            os.remove(local_path) # Remove the Gzip file

    ftp.quit()



if __name__ == "__main__":

    ftp_server1 = 'ftp.ebi.ac.uk'
    ftp_server2 = 'ftp.ncbi.nlm.nih.gov'

    gencode_remote_paths = [
        '/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz',
        '/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gff3.gz',
    ]

    refseq_remote_paths = [
        '/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz',
        '/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz',
    ]

    mane_remote_paths = [
        '/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.ensembl_genomic.gff.gz',
    ]

    refseq_assembly_reports_remote_paths = [
        '/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt',
        '/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt'
    ]

 
    script_directory = os.path.dirname(os.path.abspath(__file__))
    local_folder = os.path.join(script_directory, "data")

    if not os.path.exists(local_folder):
        os.makedirs(local_folder)
        print("Directory created successfully! Downloading files..")
        
    download_and_extract_files(ftp_server1, gencode_remote_paths, local_folder)
    download_and_extract_files(ftp_server2, refseq_remote_paths, local_folder)
    download_and_extract_files(ftp_server2, mane_remote_paths, local_folder)
    download_and_extract_files(ftp_server2, refseq_assembly_reports_remote_paths, local_folder)
    print("Files downloaded successfully!")

