import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import paths
import subprocess
from Bio.Blast import NCBIXML, NCBIWWW
import re
import pandas as pd
from Bio import SeqIO

def run_blast_online(query_fasta, db, output, evalue=1e-5,hitlist_size=50,entrez_query="(none)",matrix_name="BLOSUM62"):
    with open(query_fasta) as query_file:
        query_sequence = query_file.read()
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
                    hit_info_dict['sequence'] = hit_hsp.query
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
            # pdb.set_trace()
            sequence = hit_hsp.query
            if sequence == original_sequence:
                return True
            return False

def find_reciprocal_best_hits(original_protein_fasta, target_protein_db, output_folder):
    def get_sequence_from_fasta(fasta_file):
        sequence = ""
        with open(fasta_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                sequence = str(record.seq)
                break  # Assuming there is only one sequence in the fasta file
        return sequence

    original_protein_sequence = get_sequence_from_fasta(original_protein_fasta)

    xml_folder = os.path.join(output_folder, "blast_results")
    # if not os.path.exists(os.path.join(xml_folder, "initial_blast.xml")):
    run_blast_online(original_protein_fasta, target_protein_db, os.path.join(xml_folder, "initial_blast.xml"),
                      hitlist_size=1000,matrix_name = "BLOSUM80")   
    
    # Step 2: Parse initial BLAST results
    best_hits_by_organism = create_list_of_organisms_best_match(os.path.join(xml_folder, "initial_blast.xml"))
    
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
        with open(fasta_path, "w") as f:
            f.write(f">{organism}\n{sequence}\n")
        if check_if_the_hit_is_reciprocal(output_folder,fasta_path, original_protein_sequence, organism):
            reciprocal_hits.append(hit)
            print(f"{organism} is a reciprocal best hit")

    # Step 5: Save reciprocal best hits as a CSV file
    reciprocal_hits_df = pd.DataFrame(reciprocal_hits)
    reciprocal_hits_csv = os.path.join(output_folder, "reciprocal_best_hits.csv")
    reciprocal_hits_df.to_csv(reciprocal_hits_csv, index=False)
        


 

if __name__ == "__main__":
    # Step 1: Define the input parameters 
    juno_fasta = os.path.join(paths.juno_path, "5F4E_B.fasta")
    target_protein_db = "nr"  # Using NCBI's non-redundant protein database
    output_folder = paths.juno_orthologs_path
    orthologs = find_reciprocal_best_hits(juno_fasta, target_protein_db,output_folder)
    # create_df_of_organisms_best_match('/home/iscb/wolfson/omriyakir/Juno_Izumo/conservation_analysis/Juno/orthologs/initial_blast.xml')
    # best_hits_by_organism = create_list_of_organisms_best_match('/home/iscb/wolfson/omriyakir/Juno_Izumo/conservation_analysis/Juno/orthologs/blast_results/initial_blast.xml')