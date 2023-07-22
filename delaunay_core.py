
from scipy.spatial import Delaunay, distance
import statistics
import numpy as np
import os
import pickle

amino_acid_codes = {
    "GLY": "G",
    "ALA": "A",
    "LEU": "L",
    "MET": "M",
    "PHE": "F",
    "TRP": "W",
    "LYS": "K",
    "GLN": "Q",
    "GLU": "E",
    "SER": "S",
    "PRO": "P",
    "VAL": "V",
    "ILE": "I",
    "CYS": "C",
    "TYR": "Y",
    "HIS": "H",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "THR": "T",
    "HOH": ""
}

# Return simplices' indexes and coordinates
def build_simplices(coords):
    tri = Delaunay(coords)
    return tri.simplices   

# Build quadruplet data structure for each simplex
def build_quadruplets(simplices, amino_acids):
    quads = []
    for sim in simplices:
        quad_indices = []
        quad_residues = []
        quad_structures = []
        for index in sim:
            if amino_acids[index] in amino_acid_codes:
                quad_indices.append(index)
                quad_residues.append(amino_acid_codes[amino_acids[index]])
            #quad_structures.append(amino_acids[index]['ss'])
        quads.append(quad_residues + quad_indices) # + quad_structures)

    return quads

# Build sequential structural representation
def build_residue_strings(simplices, carbon_alphas, amino_acids,forward=True):
    strings = []
    for index in range(len(amino_acids)):
        residue_index = index
        r_simplices = [s for s in simplices if residue_index in s]
        residue_string = []
        for s in r_simplices:
            for e in s:
                if e not in residue_string:
                    if forward and e>residue_index:
                        residue_string.append(e)
                    elif not forward and e < residue_index:
                        residue_string.append(e)
        residue_string.sort()
        final_string = []

        # Process Distances
        for r in residue_string:
            distance = carbon_alphas[residue_index] - carbon_alphas[r]
            if forward:
                relative = r-residue_index
            else:
                relative = residue_index - r
            if distance <= 10:
                distance_code = 'A'
                if distance >= 5 and distance < 7.5:
                    distance_code = 'B'
                elif distance >= 7.5:
                    distance_code = 'C'
                final_string.append('{}{}'.format(str(relative),distance_code))
        final_string = ' '.join(final_string)
        strings.append(final_string)

    return '|'.join(strings)


# Calculate edge distances for each simplex
def calculate_edges(points):
    edges = []
    for simplex in points:
        simplex_edges = []
        for i in range(len(simplex)):
            first_point = simplex[i]
            if i == len(simplex)-1:
                second_point = simplex[0]
            else:
                second_point = simplex[i+1]
            p1 = first_point[0], first_point[1], first_point[2]
            p2 = second_point[0], second_point[1], second_point[2]
            dist = distance.euclidean(p1, p2)
            simplex_edges.append(dist)
        edges.append(simplex_edges)
    return edges

# Calculate tetrahedrality for given simplex
def calculate_tetrahedrality(edges):
    l_mean = statistics.mean(edges)
    total = 0
    for i in range(4):
        for j in range(4):
            if i>j:
                total += ((edges[i] - edges[j])**2) / (15 * l_mean ** 2)
    return total

# Calculate tetrahedrality for all simplices
def calculate_all_tetrahedrality(edges):
    tetra = []
    for edge_set in edges:
        t = calculate_tetrahedrality(edge_set)
        tetra.append(t)
    return tetra

# Write quadruplets as csv file
def write_simplices(quadruplets, output_dir, simplices_name):
    simplices_filename = simplices_name + ".csv"
    output_filename = os.path.join(output_dir, simplices_filename)
    with open(output_filename, "w") as f:
        np.savetxt(f, quadruplets, fmt="%s", delimiter=',')

def read_pdb_ids(filename):
    with open(filename) as f:
        pdb_ids = f.readlines()    
    return pdb_ids

def get_log_likelihood(quadruplet, likelihoods):
    residues = quadruplet[1:5]
    residues.sort()
    key = ":".join(residues)
    if key in likelihoods:
        likelihood = likelihoods[key]
    else:
        likelihood = 0
    return likelihood

