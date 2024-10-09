import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import numpy as np
import paths
from utils import create_profile_from_msa
from sequence import Sequence_logo
from biotite.sequence.align import align_optimal,SubstitutionMatrix
from biotite.sequence import ProteinSequence
from scipy.special import rel_entr

def get_the_first_sequence_of_msa(msa_file):
    with open(msa_file) as f:
        lines = f.readlines()
    return lines[1].strip()

    # Step 2: Modify each profile
def modify_profile(profile, aligned_seq):
    new_profile = []
    seq_index = 0
    for char in aligned_seq:
        if char == '-':
            comma_column = np.zeros((profile.shape[1]))
            comma_column[21] = 1
            new_profile.append(comma_column)
        else:
            new_profile.append(profile[seq_index,:])
            seq_index += 1
    new_profile = np.array(new_profile)
    return new_profile


def create_reduced_profile(profile):
    # Define the mapping of amino acids to categories
    hydro_set = set(['G','A','V','P','L','I','M','F','W'])
    polar_set = set(['S','T','C','N','Q','Y'])
    Negative_set = set(['R','H','K'])
    Positive_set = set(['D','E'])
    letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', '-']
    reduced_profile = np.zeros(profile.shape)
    
    def amino_acid_category(amino_acid):
        if amino_acid in hydro_set:
            return 'G'
        elif amino_acid in polar_set:
            return 'Q'
        elif amino_acid in Negative_set:
            return 'R'
        elif amino_acid in Positive_set:
            return 'D'
        return amino_acid
        
    for i in range(profile.shape[0]):
        for j in range(profile.shape[1]):
            amino_acid = letters[j]
            category = amino_acid_category(amino_acid)
            category_index = letters.index(category)
            reduced_profile[i,category_index] += profile[i,j]

    return reduced_profile

def create_modified_profiles():
    msa_orthologs = os.path.join(paths.juno_orthologs_path,'MSA.aln') 
    msa_paralogs = os.path.join(paths.juno_paralogs_path,'MSA.aln')
    juno_seq_with_gaps = get_the_first_sequence_of_msa(msa_orthologs)
    folate_seq_with_gaps = get_the_first_sequence_of_msa(msa_paralogs)
    profile_orthologs = create_profile_from_msa(msa_orthologs)
    non_gap_indices = np.array([i for i, char in enumerate(juno_seq_with_gaps) if char != '-'])
    sub_profile_orthologs = profile_orthologs[non_gap_indices,:]
    profile_paralogs = create_profile_from_msa(msa_paralogs)
    non_gap_indices = np.array([i for i, char in enumerate(folate_seq_with_gaps) if char != '-'])
    sub_profile_paralogs = profile_paralogs[non_gap_indices,:]
    print(f'orthologs sub_profile shape is {sub_profile_orthologs.shape}')
    print(f'paralogs sub_profile shape is {sub_profile_paralogs.shape}')
    juno_seq = ProteinSequence(juno_seq_with_gaps.replace('-',''))
    folate_seq = ProteinSequence(folate_seq_with_gaps.replace('-',''))
    alph = ProteinSequence.alphabet
    matrix = SubstitutionMatrix(alph, alph, "BLOSUM62")
    gap_penalty = -15
    alignment = align_optimal(juno_seq, folate_seq, matrix, gap_penalty=gap_penalty)[0]
    aligned_juno, aligned_folate = alignment.get_gapped_sequences()[0],alignment.get_gapped_sequences()[1]
    modified_orthologs_profile = modify_profile(sub_profile_orthologs, aligned_juno)
    print(f'modifies orthologs shape is {modified_orthologs_profile.shape}')
    modified__paralogs_profile = modify_profile(sub_profile_paralogs, aligned_folate)
    print(f'modifies paralogs shape is {modified__paralogs_profile.shape}')
    modified_orthologs_profile = create_reduced_profile(modified_orthologs_profile)
    modified__paralogs_profile = create_reduced_profile(modified__paralogs_profile)
    return modified_orthologs_profile, modified__paralogs_profile,aligned_juno, aligned_folate


def create_combined_profile_image(cavity_indexes,reduced=False):
    modified_orthologs_profile, modified__paralogs_profile,aligned_juno, aligned_folate = create_modified_profiles()
    if reduced:
        modified_orthologs_profile = create_reduced_profile(modified_orthologs_profile)
        modified__paralogs_profile = create_reduced_profile(modified__paralogs_profile)
    
    def create_real_indices_vector(aligned_seq):
        real_indices_vector = []
        cnt = 1
        for char in aligned_seq:
            if char != '-':
                real_indices_vector.append(cnt)
                cnt += 1
            else:
                real_indices_vector.append(0)
        return real_indices_vector
    real_indices_vector = [create_real_indices_vector(aligned_juno)]


    fig, (ax_indices, ax_orthologs, ax_paralogs) = plt.subplots(3, 1, figsize=(150, 11), gridspec_kw={'height_ratios': [1, 5, 5]}) 
    table = ax_indices.table(cellText=real_indices_vector, loc='center',cellLoc='center') 
    ax_indices.axis('off')
    for key, cell in table.get_celld().items():
        cell.set_text_props(weight='bold')
    for i in range(len(real_indices_vector[0])):
        if real_indices_vector[0][i] in cavity_indexes:  # Example condition: mark even numbers
            table[(0, i)].set_facecolor('yellow')

    Sequence_logo(modified_orthologs_profile, ax=ax_orthologs, ylabel="Orthologs", title="Orthologs Sequence Logo")
    Sequence_logo(modified__paralogs_profile, ax=ax_paralogs, ylabel="Paralogs", title="Paralogs Sequence Logo")

    plt.tight_layout()
    fig_path = os.path.join(paths.juno_conservation_path, f'combined_profiles_penalty_{abs(gap_penalty)}.png')
    if reduced:
        fig_path = os.path.join(paths.juno_conservation_path, f'combined_profiles_penalty_{abs(gap_penalty)}_reduced.png')
    plt.savefig(fig_path)

    
def calculate_dkl(p, q, epsilon=1e-10):
    """Calculate the Kullback-Leibler divergence between two profiles."""
    p = np.clip(p, epsilon, None)
    q = np.clip(q, epsilon, None)
    return np.sum(rel_entr(p, q))

def measure_dkl(profiles, indices):
    profile1, profile2 = profiles
    indices_set = set(indices)
    
    # Extract sub-profiles based on the given indices
    sub_profile1_in = profile1[indices, :]
    sub_profile2_in = profile2[indices, :]
    
    sub_profile1_out = np.delete(profile1, indices, axis=0)
    sub_profile2_out = np.delete(profile2, indices, axis=0)
    
    # Normalize the sub-profiles by the number of indices
    sub_profile1_in = sub_profile1_in / len(indices)
    sub_profile2_in = sub_profile2_in / len(indices)
    
    sub_profile1_out = sub_profile1_out / (profile1.shape[0] - len(indices))
    sub_profile2_out = sub_profile2_out / (profile1.shape[0] - len(indices))
    
    # Compute the Dkl for the sub-profiles
    dkl_in = calculate_dkl(sub_profile1_in, sub_profile2_in)
    dkl_out = calculate_dkl(sub_profile1_out, sub_profile2_out)
    
    return dkl_in, dkl_out


if __name__ == "__main__":   
    cavity_indexes = [28,29,30,31,32,33,64,66,67,68,69,72,73,74,75,76,77,87,90,91,93,94,97,98,110,123,124,125,126,127,139,153,154,155,157,159,188,189,190,191,192,193,194]
    # create_combined_profile_image(cavity_indexes,reduced=True)
    modified_orthologs_profile, modified__paralogs_profile,aligned_juno, aligned_folate = create_modified_profiles()
    profiles = (modified_orthologs_profile, modified__paralogs_profile)
    dkl_in, dkl_out = measure_dkl(profiles, cavity_indexes)
    print(f'DKL for the cavity indices is {dkl_in}')
    print(f'DKL for the non cavity indices is {dkl_out}')
    print(f'DKL ratio is {dkl_in/dkl_out}')
