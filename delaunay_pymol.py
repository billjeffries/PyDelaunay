from . import delaunay_core as core

# Create simplex lines between carbon alphas
def draw_lines(points, CYLINDER):
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
def select_alpha_carbons(selection, cmd):
    cmd.select("alphas", "bca.  {}".format(selection))
    atoms = cmd.get_model("alphas")
    residues = []
    ca_indices = []
    for at in atoms.atom:
        residues.append(at.resn)
        ca_indices.append(at.index)
        '''
        residues.append({"amino_acid":at.resn, 
                            "atom_index": at.index,
                            "x": at.coord[0],
                            "y": at.coord[1],
                            "z": at.coord[2],
                            "ss": "C"
                        })
        '''
    return residues, ca_indices

# Assign secondary structure to each amino acid
def assign_secondary_structures(name, ss_type, amino_acids, cmd):
    cmd.select(name, "bca. ss {}".format(ss_type))
    ss = cmd.get_model(name)
    for atom in ss.atom:
        for aa in amino_acids:
            if aa["atom_index"] == atom.index:
                aa["ss"] = ss_type 
    return amino_acids

def get_alpha_carbon_indices(atoms, ca_indices):
    indices = []
    for a in atoms:
        atom_index = a[1]
        for index, ca_index in enumerate(ca_indices):
            if ca_index == atom_index:
                indices.append(index)
    
    return indices

def get_simplices_with_amino_acids(aa_indices, simplices):
    filtered_simplices = []
    for aai in aa_indices:
        for index, s in enumerate(simplices):
            if aai in s:
                filtered_simplices.append(index)
    
    return filtered_simplices

def get_points_for_simplices(selected_simplices, points):
    filtered_points = []
    for s in selected_simplices:
        filtered_points.append(points[s])
    return filtered_points

# Return simplices' indexes and coordinates
def get_tessellation_points(selection='all', cmd=None):
    coords = cmd.get_coords(selection, 1)
    simplices = core.build_simplices(coords)
    points = coords[simplices]

    return simplices, points  


