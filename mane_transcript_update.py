import subprocess
import csv

def process_mane_gff(gff_path, shell_script_path):
    try:
        # Run the shell script and capture its output
        completed_process = subprocess.run(['bash', shell_script_path, gff_path], check=True, stdout=subprocess.PIPE, text=True)
        script_output = completed_process.stdout

        # Extracted data
        extracted_data = [line.strip().split(',') for line in script_output.split('\n')]

        # Write the extracted data to a CSV file
        csv_file_path = "./data/mane_transcripts.csv"
        with open(csv_file_path, "w", newline="") as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(["gene_name", "transcript_id", "transcript_id_refseq"])
            csv_writer.writerows(extracted_data)

    except subprocess.CalledProcessError as e:
        print(f"Script execution failed with error: {e}")
    


        




