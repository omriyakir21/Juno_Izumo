import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import paths
import subprocess
from Bio.Blast import NCBIXML, NCBIWWW
import re
import pandas as pd
from utils import align_sequences, get_sequence_from_fasta, create_profile_from_msa,save_profile_to_csv,logmaker
import biotite.sequence as seq
import biotite.application.muscle as muscle
import biotite.sequence.graphics as graphics
import matplotlib.pyplot as plt
from sequence import Sequence_logo
import numpy as np
def run_blast_online(query_sequence, db, output, evalue=1e-5,hitlist_size=50,entrez_query="(none)",matrix_name="BLOSUM62"):
    result_handle = NCBIWWW.qblast("blastp", db, query_sequence, expect=evalue, format_type="XML",
                                   hitlist_size=hitlist_size,entrez_query=entrez_query,matrix_name=matrix_name)
    with open(output, "w") as out_handle:
        out_handle.write(result_handle.read())


def run_blast_local(query_fasta, db, output, evalue=1e-5, hitlist_size=50, matrix_name="BLOSUM62"):
    # blastp_cline = NcbiblastpCommandline(query=query_fasta, db=db, evalue=evalue, outfmt=5, out=output,remote=False,
    #                                      max_target_seqs=hitlist_size, matrix=matrix_name)
    # stdout, stderr = blastp_cline()
    # if stderr:
    #     raise Exception(f"BLAST error: {stderr}")
    # with open(output, "w") as out_handle:
    #     out_handle.write(stdout)
    blast_command = [
    "blastp",
    "-query", query_fasta,       
    "-db", db,  
    "-evalue", str(evalue),              
    "-outfmt", "5", 
    "-max_target_seqs",str(hitlist_size),
    "-matrix", matrix_name,               
    "-out", output]       


    try:
        subprocess.run(blast_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("BLAST search completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e.stderr.decode()}")


def create_list_of_organisms_best_match(blast_xml):
    import pdb 
    # pdb.set_trace()
    organisms_set = set()
    with open(blast_xml) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        best_hits_by_organism = []
        # pdb.set_trace()
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                hit_description = alignment.hit_def
                 # Use regex to extract organism name (assuming it is in brackets)
                match = re.search(r'\[(.*?)\]', hit_description)
                if match:
                    organism_name = match.group(1)
                    if organism_name in organisms_set:
                        continue
                    organisms_set.add(organism_name)
                    hit_info_dict = {}
                    hit_info_dict['organism'] = match.group(1)
                    hit_hsp = alignment.hsps[0]
                    hit_info_dict['evalue'] = hit_hsp.expect
                    hit_info_dict['sequence'] = hit_hsp.sbjct
                    hit_info_dict['title'] = hit_description
                    best_hits_by_organism.append(hit_info_dict)
   
    return best_hits_by_organism
    

def check_if_the_hit_is_reciprocal(output_folder,query_fasta,original_sequence,organism):
    # import pdb
    xml_folder = os.path.join(output_folder, "blast_results")
    reciprocal_blast_xml = os.path.join(xml_folder,f"{organism}_reciprocal_blast.xml")
    run_blast_local(query_fasta, paths.human_proteins_db_path, reciprocal_blast_xml,matrix_name="BLOSUM80",hitlist_size=5)
    with open(reciprocal_blast_xml) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            alignment = blast_record.alignments[0]
            hit_hsp = alignment.hsps[0]
            sequence = hit_hsp.sbjct
            if sequence == original_sequence or align_sequences(sequence,original_sequence) > 0.95:
                return True
            return False

def run_inital_blast(original_protein_fasta, db, output_folder,hitlist_size,matrix_name = "BLOSUM80"):
    original_protein_sequence = get_sequence_from_fasta(original_protein_fasta)
    xml_folder = os.path.join(output_folder, "blast_results")
    if not os.path.exists(xml_folder):
        os.mkdir(xml_folder)
    run_blast_online(original_protein_sequence, db, os.path.join(xml_folder, f"initial_blast_{hitlist_size}.xml"),hitlist_size = hitlist_size,matrix_name = matrix_name)

def find_reciprocal_best_hits(original_protein_fasta, target_protein_db, output_folder,hitlist_size):
    original_protein_sequence = get_sequence_from_fasta(original_protein_fasta)
    xml_folder = os.path.join(output_folder, "blast_results")

    # Step 2: Parse initial BLAST results
    best_hits_by_organism = create_list_of_organisms_best_match(os.path.join(xml_folder, f"initial_blast_{hitlist_size}.xml"))
    
    # Step 3: Save the DataFrame as a CSV file
    best_matches_df = pd.DataFrame(best_hits_by_organism)
    output_csv = os.path.join(output_folder, "best_matches.csv")
    best_matches_df.to_csv(output_csv, index=False)

    fasta_files_to_align_dir = os.path.join(output_folder, 'fasta_files_to_align')
    if not os.path.exists(fasta_files_to_align_dir):
        os.mkdir(fasta_files_to_align_dir)
    # Step 4: Filter for reciprocal best hits
    reciprocal_hits = []
    for hit in best_hits_by_organism:
        organism = hit['organism']
        sequence = hit['sequence']

        fasta_path = os.path.join(fasta_files_to_align_dir, f"{organism}.fasta")
        if not os.path.exists(fasta_path):
            with open(fasta_path, "w") as f:
                f.write(f">{organism}\n{sequence}\n")
        if check_if_the_hit_is_reciprocal(output_folder,fasta_path, original_protein_sequence, organism):
            reciprocal_hits.append(hit)
            print(f"{organism} is a reciprocal best hit")

    # Step 5: Save reciprocal best hits as a CSV file
    reciprocal_hits_df = pd.DataFrame(reciprocal_hits)
    reciprocal_hits_csv = os.path.join(output_folder, "reciprocal_best_hits.csv")
    reciprocal_hits_df.to_csv(reciprocal_hits_csv, index=False)
        

def get_sequences_from_csv(csv_path):
    df = pd.read_csv(csv_path)
    sequences = df['sequence'].tolist()
    sequences_without_gaps = [re.sub(r'[^A-Z]', '', sequence) for sequence in sequences]
    return sequences_without_gaps
 
def get_organisms_from_csv(csv_path):
    df = pd.read_csv(csv_path)
    organisms = df['organism'].tolist()
    return organisms

def construct_MSA(sequences,organisms, output_folder):
    assert len(sequences) == len(organisms)
    sequences = [seq.ProteinSequence(sequence) for sequence in sequences]
    app = muscle.Muscle5App(sequences)
    app.start()
    app.join()
    alignment = app.get_alignment()
    gapped_seqs = alignment.get_gapped_sequences()
    for i in range(len(organisms)):
        print(organisms[i]," "*3,gapped_seqs[i])
        # Save the alignment to a file in FASTA format
    with open(os.path.join(output_folder, "orthologs_alignment_muscle5.fasta"), "w") as output_file:
        for i, gapped_seq in enumerate(gapped_seqs):
            output_file.write(f">{organisms[i]}\n{gapped_seq}\n")
    fig = plt.figure(figsize=(8.0, 8.0))
    ax = fig.add_subplot(111)
    order = app.get_alignment_order()
    graphics.plot_alignment_type_based(
    ax,
    alignment[:200, order.tolist()],labels=[organisms[i] for i in order],show_numbers=True,color_scheme="clustalx")
    fig.tight_layout()
    plt.savefig(os.path.join(output_folder, "MSA.png"))
    output_path = os.path.join(output_folder, "MSA.aln")
    with open(output_path, "w") as f:
        for i in range(len(organisms)):
            f.write(f">{organisms[i]}\n{gapped_seqs[i]}\n")

    return output_path




if __name__ == "__main__":
    # Step 1: Define the input parameters 
    # juno_fasta = os.path.join(paths.juno_data_path, "5F4E_B.fasta")
    izumo_fasta = os.path.join(paths.izumo_data_path, "5F4E_A.fasta")
    target_protein_db = "nr"  # Using NCBI's non-redundant protein database
    # output_folder = paths.juno_orthologs_path
    output_folder = paths.izumo_orthologs_path
    hitlist_size = 5000
    run_inital_blast = run_inital_blast(izumo_fasta, target_protein_db, output_folder,hitlist_size)
    # orthologs = find_reciprocal_best_hits(izumo_fasta, target_protein_db,output_folder,hitlist_size)
    # create_df_of_organisms_best_match('/home/iscb/wolfson/omriyakir/Juno_Izumo/conservation_analysis/Juno/orthologs/initial_blast.xml')
    # best_hits_by_organism = create_list_of_organisms_best_match('/home/iscb/wolfson/omriyakir/Juno_Izumo/conservation_analysis/Juno/orthologs/blast_results/initial_blast.xml')
    # sequences = get_sequences_from_csv(os.path.join(output_folder, "reciprocal_best_hits.csv"))
    # organisms = get_organisms_from_csv(os.path.join(output_folder, "reciprocal_best_hits.csv"))
    # msa = construct_MSA(sequences,organisms,output_folder)
    # profile = create_profile_from_msa(os.path.join(output_folder,'MSA.aln'))
    # with open(os.path.join(output_folder,'MSA.aln'), 'r') as msa_file:
    #     first_row = msa_file.readline().strip()
    #     second_row = msa_file.readline().strip()
    #     non_gap_indices = np.array([i for i, char in enumerate(second_row) if char != '-'])
    
    # sub_profile = profile[non_gap_indices,:]
    # # Step 6: Save the profile as a sequence logo
    # sequence_logo = Sequence_logo(sub_profile)
    # sequence_logo.savefig(os.path.join(output_folder, 'sequence_logo_sub_profile.png'))


