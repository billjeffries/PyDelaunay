from pymol import cgo, cmd
from pymol.cgo import *    

from scipy.spatial import Delaunay, distance
import statistics
import numpy as np
import os

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
    "THR": "T"
}

# Return simplices' indexes and coordinates
def get_tessellation_points(selection='all'):
    coords = cmd.get_coords(selection, 1)
    tri = Delaunay(coords)
    tri_points = coords[tri.simplices]
    return tri.simplices, tri_points  

# Create simplex lines between carbon alphas
def draw_lines(points):
    first_point = points[0][0]
    second_point = points[0][1]
    obj = [
        CYLINDER,  first_point[0], first_point[1], first_point[2], 
            second_point[0], second_point[1], second_point[2], 
            0.05,
            1.0, 0.0, 0.0,
            1.0, 0.0,0.0,
    ]
    for simplex in points:
        for i in range(len(simplex)):
            if i == len(simplex)-1:
                break
            first_point = simplex[i]
            second_point = simplex[i+1]
            obj.extend([
                CYLINDER,  first_point[0], first_point[1], first_point[2], 
                    second_point[0], second_point[1], second_point[2], 
                    0.05,
                    1.0, 0.0, 0.0,
                    1.0, 0.0,0.0,
            ])
    return obj

# Select carbon alphas and return info as array of carbon alpha dictionaries
def select_alpha_carbons(selection):
    cmd.select("alphas", "bca.  {}".format(selection))
    atoms = cmd.get_model("alphas")
    residues = []
    for at in atoms.atom:
        residues.append({"amino_acid":at.resn, 
                            "atom_index": at.index,
                            "x": at.coord[0],
                            "y": at.coord[1],
                            "z": at.coord[2],
                            "ss": "C"
                        })
        '''
        print("ATOM DEFINITION: "+at.chain+" "\
                                +at.resn+" "\
                                +str(at.resi)+" "\
                                +at.name+" "\
                                +str(at.index)+" "\
                                +"%.2f " % (at.b)\
                                +str(at.coord[0])+" "\
                                +str(at.coord[1])+" "\
                                +str(at.coord[2]))
        '''
    return residues

# Assign secondary structure to each amino acid
def assign_secondary_structures(name, ss_type, amino_acids):
    cmd.select(name, "bca. ss {}".format(ss_type))
    ss = cmd.get_model(name)
    for atom in ss.atom:
        for aa in amino_acids:
            if aa["atom_index"] == atom.index:
                aa["ss"] = ss_type 
    return amino_acids

# Build quadruplet data structure for each simplex
def build_quadruplets(simplices, amino_acids):
    quads = []
    for sim in simplices:
        quad_indices = []
        quad_residues = []
        quad_structures = []
        for index in sim:
            quad_indices.append(index)
            quad_residues.append(amino_acid_codes[amino_acids[index]['amino_acid']])
            quad_structures.append(amino_acids[index]['ss'])
        quads.append(quad_residues + quad_indices + quad_structures)

    return quads

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
def write_simplices(quadruplets, output_file):
    with open(output_file, "a") as f:
        #numpy.savetxt(f, a)    
        np.savetxt(f, quadruplets, fmt="%s", delimiter=',')

def get_alpha_carbon_indices(atoms):
    indices = []
    for a in atoms:
        atom_index = a[1]
        for index, aa in enumerate(amino_acids):
            if aa['atom_index'] == atom_index:
                indices.append(index)
    
    return indices

def get_simplices_with_amino_acids(aa_indices):
    filtered_simplices = []
    for aai in aa_indices:
        for index, s in enumerate(simplices):
            if aai in s:
                filtered_simplices.append(index)
    
    return filtered_simplices

def get_points_for_simplices(selected_simplices):
    filtered_points = []
    for s in selected_simplices:
        filtered_points.append(points[s])
    return filtered_points

def read_pdb_ids(filename):
    with open(filename) as f:
        pdb_ids = f.readlines()    
    return pdb_ids

# PyMOL commands
@cmd.extend
def delaunay(selection='all', output_file=None, name='delaunay'):
    global simplices, points, amino_acids
    amino_acids = select_alpha_carbons(selection)
    simplices, points = get_tessellation_points("alphas")
    lines = draw_lines(points)
    amino_acids = assign_secondary_structures("helices", "H", amino_acids)
    amino_acids = assign_secondary_structures("sheets", "S", amino_acids)
    quadruplets = build_quadruplets(simplices, amino_acids)
    edges = calculate_edges(points)
    tetrahedrality = calculate_all_tetrahedrality(edges)
    for i in range(len(quadruplets)):
        quadruplets[i].insert(0, selection)
        quadruplets[i].append(tetrahedrality[i])
    if output_file is not None:
        write_simplices(quadruplets, output_file)

    cmd.delete("alphas")
    cmd.delete("helices")
    cmd.delete("sheets")
    cmd.load_cgo(lines,name) 

@cmd.extend
def delaunay_filter(selection='all', name='delaunay_filter'):
    atoms = cmd.index(selection)
    selected_carbon_alphas = get_alpha_carbon_indices(atoms)
    filtered_simplices = get_simplices_with_amino_acids(selected_carbon_alphas)
    filtered_points = get_points_for_simplices(filtered_simplices)
    filtered_lines = draw_lines(filtered_points)
    cmd.disable('delaunay')
    cmd.load_cgo(filtered_lines,name) 

@cmd.extend
def delaunay_batch(input_filename, output_filename):
    pdb_ids = read_pdb_ids(input_filename)
    for pdb_id in pdb_ids:
        sel = cmd.fetch(pdb_id, '', 0, 1, -1, -2, -1, 'pdb', 0)
        delaunay(sel, output_filename)
        cmd.delete("delaunay")
        cmd.delete(pdb_id)

    for pdb_id in pdb_ids:
        os.remove("{}.pdb".format(pdb_id.replace("\n","")))

# Map PyMOL command to function
cmd.extend('delaunay', delaunay)
cmd.extend('delaunay_filter', delaunay_filter)
cmd.extend('delaunay_batch', delaunay_batch)