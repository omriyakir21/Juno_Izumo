import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import subprocess
from Bio.Blast import NCBIXML
import paths
from conservation_analysis.find_orthologs import get_sequences_from_csv,get_organisms_from_csv,run_blast_online,construct_MSA
import math
from Bio import pairwise2
from utils import cluster_sequences
import logging
from utils import align_sequences, get_sequence_from_fasta, create_profile_from_msa,save_profile_to_csv,logmaker
import numpy as np
import pandas as pd
from sequence import Sequence_logo


def find_paralogs_for_orthologs(sequences, organisms,output_folder):
    for sequence, organism in zip(sequences, organisms):
        paralogs = find_paralogs(sequence, organism, output_folder)
        if paralogs is None:
            logging.info(f"Paralogs for {organism} already exist")
            continue
        logging.info(f"Paralogs for {organism}:")
        for paralog in paralogs:
            logging.info(f"  {paralog['aligned_sequence']} ({paralog['similarity']:.0%})")
        logging.info("")
    print("Done")


def is_unique_paralog(results, gapped_sequence):
    sequence = gapped_sequence.replace("-", "")
    for result in results:
        aligned_seq = result['aligned_sequence'].replace("-", "")
        alignments = pairwise2.align.globalxx(sequence, aligned_seq)
        max_identity = max(
            sum(1 for a, b in zip(aln[0], aln[1]) if a == b) / len(aln[0])
            for aln in alignments
        )
        if max_identity >= 0.95:
            return False
    return True


def find_paralogs(sequence, organism, output_folder):
    output_file = os.path.join(output_folder, f"{organism}_paralogs.csv")
    #check if outputfile exists and not empty
    # if os.path.exists(output_file):
    #     return None
    # Run BLAST command with entrez_query to filter by organism
    output_path = os.path.join(paths.izumo_paralogs_blast_results_path, f'{organism}.xml')
    if not os.path.exists(output_path):
        attempts = 0
        max_attempts = 3
        while attempts < max_attempts:
            try:
                run_blast_online(sequence, "nr", output_path, entrez_query=f'{organism}[Organism] OR {organism}[ORGN]')
                break
            except Exception as e:
                attempts += 1
                if attempts == max_attempts:
                    logging.error(f"Failed to run BLAST online after {max_attempts} attempts, sequence = {sequence}, organism = {organism}")
                    return None


    # Parse BLAST results
    with open(output_path) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        
        all_sequences = []
        
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    similarity = hsp.identities / hsp.align_length
                    all_sequences.append({
                        "aligned_sequence": hsp.sbjct,
                        "definition": alignment.hit_def,
                        "similarity": similarity
                    })

        paralogs = []
        if all_sequences:
            # Cluster the sequences and take cluster representatives
            _, representative_indices = cluster_sequences([all_sequences[i]["aligned_sequence"].replace("-", "") for i in range(len(all_sequences))],
                                                            seqid=0.95,coverage=0.9, covmode='0')
            
            paralogs = [all_sequences[i] for i in representative_indices if all_sequences[i]['similarity'] <= 0.95]
            


        # Save the cluster representatives as the final paralogs
        
        with open(output_file, "w") as csvfile:
            csvfile.write("sequence,similarity,definition\n")
            for result in paralogs:
                csvfile.write(f"{result['aligned_sequence']},{result['similarity']:.2f},{result['definition']}\n")

    return paralogs

import os
import pandas as pd

def process_paralogs_folder(folder_path, output_folder):
    msa = os.path.join(output_folder,'MSA.aln')
    if not os.path.exists(msa):

        # Step 1: Read and concatenate all CSV files in the folder
        all_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.csv')]
        
        df_list = []
        for file in all_files:
            try:
                print(f"Reading {file}")
                with open(file, 'r') as f:
                    lines = f.readlines()
                    
                # Process each line to replace redundant commas and remove '-' from sequences
                processed_lines = []
                for line in lines:
                    parts = line.strip().split(',')
                    if len(parts) > 3:
                        parts[2] = '_'.join(parts[2:])
                        parts = parts[:3]
                    parts[0] = parts[0].replace('-', '')
                    processed_lines.append(','.join(parts))
                
                # Write the processed lines to a temporary file and read it with pandas
                temp_file = os.path.join(output_folder, 'temp.csv')
                with open(temp_file, 'w') as f:
                    f.write('\n'.join(processed_lines))
                
                df = pd.read_csv(temp_file)
                
                # Add the organism column
                organism_name = os.path.basename(file).replace('_paralogs.csv', '')
                df['organism'] = organism_name
                
                df_list.append(df)
            except Exception as e:
                print(f"Failed to read {file}: {e}")
                raise
        
        concatenated_df = pd.concat(df_list, ignore_index=True)
        
        # Step 2: Save concatenated data to a new CSV file
        output_csv_path = os.path.join(output_folder, 'Homo sapiens_paralogs.csv')
        concatenated_df.to_csv(output_csv_path, index=False)

        # Step 3: Extract sequences and organisms from the concatenated data
        sequences = concatenated_df['sequence'].tolist()
        organisms = concatenated_df['organism'].tolist()
        # Step 4: Construct MSA
        msa = construct_MSA(sequences,organisms,output_folder)
    
    # Step 5: Create a profile from the MSA
    profile = create_profile_from_msa(msa)
    # Create a sub-profile containing only the columns of the non '-' letters of the first row of the MSA
    with open(msa, 'r') as msa_file:
        first_row = msa_file.readline().strip()
        second_row = msa_file.readline().strip()
        non_gap_indices = np.array([i for i, char in enumerate(second_row) if char != '-'])
    
    sub_profile = profile[non_gap_indices,:]
    # Step 6: Save the profile as a sequence logo
    sequence_logo = Sequence_logo(sub_profile)
    sequence_logo.savefig(os.path.join(output_folder, 'sequence_logo_sub_profile.png'))




if __name__ == "__main__":
    input_folder = paths.izumo_orthologs_path
    # input_folder = paths.juno_orthologs_path
        # Get the batch index from command line arguments
    sequences = get_sequences_from_csv(os.path.join(input_folder, "reciprocal_best_hits.csv"))
    organisms = get_organisms_from_csv(os.path.join(input_folder, "reciprocal_best_hits.csv"))
    # Configure logging

    # # batch_index = int(sys.argv[1])
    # batch_index = 2
    # batch_size = 200
    # start_index = batch_index * batch_size
    # end_index = min(start_index + batch_size, len(sequences))
    # sequences = sequences[start_index:end_index]
    # organisms = organisms[start_index:end_index]
    # output_folder = paths.juno_organisms_paralogs_csvs_path
    output_folder = paths.izumo_organisms_paralogs_csvs_path
    log_file = os.path.join(paths.izumo_paralogs_path, f"find_paralogs.log")
    logging.basicConfig(filename=log_file, level=logging.INFO, 
                        format='%(asctime)s %(levelname)s:%(message)s')
    # logging.info("Starting find_paralogs script", extra={'batch_index': batch_index})
    find_paralogs_for_orthologs(sequences, organisms,output_folder)
    # process_paralogs_folder(output_folder,paths.juno_paralogs_path)

