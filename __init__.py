import warnings
from . import delaunay_core as core, delaunay_pymol as dpymol

try:
    from pymol import cgo, cmd
    from pymol.cgo import *    
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: Ignore if you are only using pdb reader module')
    cmd = None

def __init_plugin__(app=None):
    if cmd is not None:
        pass

if cmd is not None:
    # PyMOL commands
    @cmd.extend
    def delaunay(selection='all'):
        global simplices, points, amino_acids, ca_indices

        cmd.delete("delaunay")

        amino_acids, ca_indices = dpymol.select_alpha_carbons(selection, cmd)
        simplices, points = dpymol.get_tessellation_points("alphas", cmd)
        lines = dpymol.draw_lines(points, CYLINDER)
        quadruplets = core.build_quadruplets(simplices, amino_acids)

        output_file = '/Users/billjeffries/'.format(selection)
        core.write_simplices(quadruplets, output_file, selection)

        cmd.delete("alphas")
        cmd.delete("helices")
        cmd.delete("sheets")
        cmd.load_cgo(lines,'delaunay') 

    @cmd.extend 
    def delaunay_export(selection="all", output_dir=None):
        amino_acids, ca_indices = dpymol.select_alpha_carbons(selection, cmd)
        simplices, _ = dpymol.get_tessellation_points("alphas", cmd)
        quadruplets = core.build_quadruplets(simplices, amino_acids)

        output_file = '{}/{}.csv'.format(output_dir, selection)
        core.write_simplices(quadruplets, output_dir, selection)
        print("{} exported to {}".format(selection, output_file))

        cmd.delete("alphas")

    @cmd.extend
    def delaunay_export_by_residue(selection="all", output_dir=None):
        amino_acids, ca_indices = dpymol.select_alpha_carbons(selection, cmd)
        simplices, _ = dpymol.get_tessellation_points("alphas", cmd)

    @cmd.extend
    def delaunay_filter(selection='all', name='delaunay_filter'):
        cmd.delete("delaunay_filter")

        atoms = cmd.index(selection)
        selected_carbon_alphas = dpymol.get_alpha_carbon_indices(atoms, ca_indices)
        filtered_simplices = dpymol.get_simplices_with_amino_acids(selected_carbon_alphas,simplices)
        filtered_points = dpymol.get_points_for_simplices(filtered_simplices, points)
        filtered_lines = dpymol.draw_lines(filtered_points, CYLINDER)
        cmd.disable('delaunay')
        cmd.load_cgo(filtered_lines,name) 

