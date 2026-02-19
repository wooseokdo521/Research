import os
import subprocess
import pandas as pd
from Bio import SeqIO

#path
reference_seq = "/home/xtremodevolab/Desktop/Daniel/medaka_HEs.fasta"
genome_folder = "/home/xtremodevolab/Desktop/Daniel/Genome"
output_folder = "/home/xtremodevolab/Desktop/Daniel/BLAST"

os.makedirs(output_folder, exist_ok=True)

#loop all genomes for blast
for genome_file in os.listdir(genome_folder):
    if genome_file.endswith(".faa"):
        genome_path = os.path.join(genome_folder, genome_file)
        
        #make blast db
        subprocess.run([
            "makeblastdb",
            "-in", genome_path,
            "-dbtype", "prot"
        ])
        
        #set output
        if "_protein.faa" in genome_file:
            base_name = genome_file.replace("_protein.faa", "")
        else:
            base_name = genome_file.rsplit(".", 1)[0]  # removes .faa
        
        raw_output_file = os.path.join(output_folder, base_name + "_blastp.tsv")
        filtered_output_file = os.path.join(output_folder, base_name + "_blastp_filtered.tsv")

        #Run blastp
        subprocess.run([
            "blastp",
            "-query", reference_seq,
            "-db", genome_path,
            "-out", raw_output_file,
            "-outfmt", "6",
            "-evalue", "1e-5",
            "-num_threads", "4"
        ])
        
        #headers
        headers = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
        data = pd.read_csv(raw_output_file, sep="\t", header=None, names=headers)

        data.to_csv(raw_output_file, sep="\t", index=False)

        #filter top 100 hits
        data_filtered = data.sort_values("bitscore", ascending=False).head(100)
        data_filtered.to_csv(filtered_output_file, sep="\t", index=False)



        print(f"BLAST completed for {genome_file}")

print("All BLAST searches completed.")

#parse all top 100 hits
dfs = []

#loop all filtered files
for file in os.listdir(output_folder):
    if file.endswith("_blastp_filtered.tsv"):
        file_path = os.path.join(output_folder, file)

        df = pd.read_csv(file_path, sep="\t")

        #add column for genome name
        genome_name = file.replace("_blastp_filtered.tsv", "")
        df["genome"] = genome_name

        dfs.append(df)

#combine all into one file
combined_df = pd.concat(dfs, ignore_index=True)
selected_cols = ["genome","sseqid","pident","length","sstart","send","evalue"]
combined_df = combined_df[selected_cols]

parsed_file = os.path.join(output_folder, "combined_filtered_parsed_hits.tsv")
combined_df.to_csv(parsed_file, sep="\t", index=False)

print("All top 100 hits combined.")

#set output for mafft
mafft_input_data = os.path.join(output_folder, "combined_filtered_parsed_hits.faa")
mafft_output_data = os.path.join(output_folder, "combined_top100_hits_alignment.faa")

#load all sequences into dictionary
genome_seqs = {}
for genome_file in os.listdir(genome_folder):
    if genome_file.endswith(".faa"):
        genome_path = os.path.join(genome_folder, genome_file)
        for record in SeqIO.parse(genome_path, "fasta"):
            genome_seqs[record.id.split()[0]] = record

#retrieve hit sequences to fasta
with open(mafft_input_data, "w") as f:
    for sseqid in combined_df["sseqid"]:
        if sseqid in genome_seqs:
            SeqIO.write(genome_seqs[sseqid], f, "fasta")
        else:
            print(f"{sseqid} not found.")

#run mafft
subprocess.run([
    "mafft",
    "--localpair", #L-INS-i algorithm
    "--maxiterate", "1000",
    "--op", "1", #gap open penalty
    "--ep", "0.1", #offset for gap extension
    mafft_input_data
], stdout=open(mafft_output_data, "w"))

print("MAFFT alignment completed.")