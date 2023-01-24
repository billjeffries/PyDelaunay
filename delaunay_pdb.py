import os
import sys
import urllib.request
import Bio
import Bio.PDB
import Bio.SeqRecord
import delaunay_core as core
import urllib.request, json 


def process_pdb(pdb_id, download_dir):
    pdb_filename = download_pdb(pdb_id, download_dir)
    if pdb_filename is None:
        return None
    struct = read_pdb(pdb_id, pdb_filename)
    os.remove(pdb_filename)
    return struct

def download_pdb(pdb_id, download_dir, download_url="http://files.rcsb.org/download/"):
    pdb_filename = pdb_id + ".pdb"
    url = download_url + pdb_filename
    output_filename = os.path.join(download_dir, pdb_filename)
    try:
        urllib.request.urlretrieve(url, output_filename)
        return output_filename
    except Exception as err:
        # all sorts of things could have gone wrong...
        print(str(err), file=sys.stderr)
        return None

def read_pdb(pdb_id, pdb_filename):
    try:
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)  
        struct = pdb_parser.get_structure(pdb_id, pdb_filename)
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None 

def get_calphas(struct):
    calphas = [ atom for atom in struct.get_atoms() if atom.get_fullname() == " CA " ]
    return calphas

def simplices_from_pdb(pdb_id, download_dir):
    struct = process_pdb(pdb_id, download_dir)
    amino_acids = []
    for r in struct.get_residues():
        amino_acids.append(r.get_resname())

    alphas = get_calphas(struct)
    points = []
    for a in alphas:
        coords = a.coord
        points.append(coords)

    simplices = core.build_simplices(points)
    quads = core.build_quadruplets(simplices, amino_acids)
    return quads

def process_batch_pdb_simplices(pdb_ids, download_dir):
    for pdb_id in pdb_ids:
        simplices = simplices_from_pdb(pdb_id, download_dir)
        core.write_simplices(simplices, download_dir, pdb_id)

def process_all_current_pdbs(download_dir):
    with urllib.request.urlopen("https://data.rcsb.org/rest/v1/holdings/current/entry_ids") as url:
        pdb_ids = json.load(url)
        process_batch_pdb_simplices(pdb_ids, download_dir)

def process_pdbs_by_similarity(similarity, download_dir):
    pdb_filename = "pdb-cluster-" + str(similarity) + ".txt"
    output_filename = os.path.join(download_dir, pdb_filename)
    url = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-" + str(similarity) + ".txt"
    with urllib.request.urlopen(url) as f:
        clusters = f.readlines()
        for cluster in clusters:
            cluster_ids = cluster.decode('ascii')
            pdb_ids = cluster_ids.split(' ')
            print(pdb_ids)

