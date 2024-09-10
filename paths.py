import os

# List of directory and files paths
current_dir = os.path.dirname(os.path.abspath(__file__))

#data
data_path = os.path.join(current_dir, 'data')
#   Juno
juno_path = os.path.join(data_path, 'Juno')
#   Izumo
izumo_path = os.path.join(data_path, 'Izumo')
#   human_proteins_db
human_proteins_db_path = os.path.join(data_path, 'human_proteins_db')


#conservation_analysis
conservation_analysis_path = os.path.join(current_dir, 'conservation_analysis')
#   Juno
juno_conservation_path = os.path.join(conservation_analysis_path, 'Juno')
#       orthologs
juno_orthologs_path = os.path.join(juno_conservation_path, 'orthologs')
#           fasta_files_to_align
juno_fasta_files_to_align_path = os.path.join(juno_orthologs_path, 'fasta_files_to_align')
#           blast_results
juno_blast_results_path = os.path.join(juno_orthologs_path, 'blast_results')
#   Izumo
izumo_conservation_path = os.path.join(conservation_analysis_path, 'Izumo')
#       orthologs
izumo_orthologs_path = os.path.join(izumo_conservation_path, 'orthologs')
#           fasta_files_to_align
izumo_fasta_files_to_align_path = os.path.join(izumo_orthologs_path, 'fasta_files_to_align')
#           blast_results
izumo_blast_results_path = os.path.join(izumo_orthologs_path, 'blast_results')
