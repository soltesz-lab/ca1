## Instantiates a hoc cell and exports its 3d points to SWC format

import os, sys, gc, logging, random, copy, pprint
import click
import numpy as np
from mpi4py import MPI
import ca1
import ca1.io_utils as io_utils
import ca1.utils as utils
import ca1.neuron_utils as neuron_utils
from ca1.env import Env
from neuron import h, gui
from collections import defaultdict
from neuroh5.io import append_cell_attributes, append_cell_trees, read_population_ranges
import h5py
import networkx as nx

sys_excepthook = sys.excepthook
def mpi_excepthook(type, value, traceback):
    sys_excepthook(type, value, traceback)
    if MPI.COMM_WORLD.size > 1:
        MPI.COMM_WORLD.Abort(1)
sys.excepthook = mpi_excepthook


def export_swc_dict(cell, ref_axis=1, sections=[("soma_list",1),("apical_list",4),("basal_list",3),("axon_list",2)]):

    min_sec_pts = 3
    swc_point_idx = 0
    swc_points = []
    swc_point_sec_dict = defaultdict(list)
    sec_dict = {}

    origin = getattr(cell, 'soma')[0]
    if ref_axis == 0:
        origin_ref = origin.x3d(0)
    elif ref_axis == 1:
        origin_ref = origin.y3d(0)
    else:
        origin_ref = origin.z3d(0)
    
    for section, sectype in sections:
        if hasattr(cell, section):
            seclist = list(getattr(cell, section))
            for sec in seclist:
                if hasattr(sec, 'sec'):
                    sec = sec.sec
                L = sec.L
                pch_x, pch_y, pch_z, pch_diam = neuron_utils.interplocs(sec)
                npts_interp = max(int(round(L)), min_sec_pts)
                interp_locs = np.linspace(0, 1, npts_interp)
                xyz_interp = np.column_stack((pch_x(interp_locs), pch_y(interp_locs), pch_z(interp_locs), pch_diam(interp_locs)))
                sec.pt3dclear()
                h.pt3dadd(h.Vector(xyz_interp[:,0]), h.Vector(xyz_interp[:,1]),
                          h.Vector(xyz_interp[:,2]), h.Vector(xyz_interp[:,3]),
                          sec=sec)
                n3d = sec.n3d()
                assert(n3d >= min_sec_pts)
                for i in range(n3d):

                    x = sec.x3d(i)
                    y = sec.y3d(i)
                    z = sec.z3d(i)
                    d = sec.diam3d(i)
                    ll = sec.arc3d(i)
                    rad = d / 2.
                    loc = ll / L
                    distance = h.distance(origin(0.5), sec(loc))
                    coord = [x, y, z]
                    if origin_ref > coord[ref_axis]:
                        distance = -distance
                    first = True if i == 0 else False
                    layer = int(cell.get_layer(distance))
                    swc_point = (swc_point_idx, sectype, x, y, z, rad, loc, sec, first, layer)
                    swc_points.append(swc_point)
                    swc_point_sec_dict[sec.name()].append(swc_point)
                    swc_point_idx += 1

    pt_xs = []
    pt_ys = []
    pt_zs = []
    pt_radius = []
    pt_parents = []
    pt_layers = []
    pt_swc_types = []

    sec_idx_dict = {}
    sec_name_dict = {}
    for sec_idx, sec_name in enumerate(sorted(swc_point_sec_dict.keys())):
        sec_idx_dict[sec_name] = sec_idx
        sec_name_dict[sec_idx] = sec_name

    sec_graph = nx.DiGraph()
    sec_graph.add_nodes_from(range(len(sec_name_dict)))
    
    for swc_point in swc_points:
        (swc_point_idx, sectype, x, y, z, rad, loc, sec, first, layer) = swc_point
        parent_idx = -1
        if not first:
            parent_idx = swc_point_idx-1
        else:
            parent_seg = sec.parentseg()
            if parent_seg is not None:
                parent_x = parent_seg.x
                parent_sec = parent_seg.sec
                parent_points = swc_point_sec_dict[parent_sec.name()]
                parent_sec_idx = sec_idx_dict[parent_sec.name()]
                this_sec_idx = sec_idx_dict[sec.name()]
                sec_graph.add_edge(parent_sec_idx, this_sec_idx)
                for parent_point in parent_points:
                    (parent_point_idx, _, _, _, _, _, parent_point_loc, _, _, _) = parent_point
                    if parent_point_loc >= parent_x:
                        parent_idx = parent_point_idx
                        break
                parent_idx = parent_point_idx
                if parent_idx < 0:
                    raise RuntimeError(f'export_neuroh5_tree: unable to locate parent point of point {swc_point}; '
                                       f'parent segment is {parent_seg}: {pprint.pformat(parent_points)}')
        pt_xs.append(x)
        pt_ys.append(y)
        pt_zs.append(z)
        pt_radius.append(rad)
        pt_parents.append(parent_idx)
        pt_layers.append(layer)
        pt_swc_types.append(sectype)
        
    sec_src = []
    sec_dst = []

    for (src,dst) in sec_graph.edges:
        sec_src.append(src)
        sec_dst.append(dst)

    pt_sections = []
    pt_sections.append(len(sec_name_dict))
    for sec_idx in sorted(sec_name_dict.keys()):
        sec_name = sec_name_dict[sec_idx]
        sec_pts = swc_point_sec_dict[sec_name]
        pt_sections.append(len(sec_pts))
        for (swc_point_idx, sectype, x, y, z, rad, loc, sec, first, layer) in sec_pts:
            pt_sections.append(swc_point_idx)

        
    tree_dict = { 'x': np.asarray(pt_xs, dtype=np.float32),
                  'y': np.asarray(pt_ys, dtype=np.float32),
                  'z': np.asarray(pt_zs, dtype=np.float32),
                  'radius': np.asarray(pt_radius, dtype=np.float32),
                  'layer': np.asarray(pt_layers, dtype=np.int8),
                  'parent': np.asarray(pt_parents, dtype=np.int32),
                  'swc_type': np.asarray(pt_swc_types, dtype=np.int8),
                  'sections': np.asarray(pt_sections, dtype=np.uint16),
                  'src': np.asarray(sec_src, dtype=np.uint16),
                  'dst': np.asarray(sec_dst, dtype=np.uint16) }

    return tree_dict


    
@click.command()
@click.option("--config", '-c', required=True, type=str)
@click.option("--config-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), default='config')
@click.option("--population", '-p', required=True, type=str)
@click.option("--gid", '-g', default=0, type=int)
@click.option("--ref-axis", '-r', default=1, type=int)
@click.option("--template-name", '-t', required=True, type=str)
@click.option("--input-file", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--output-file", required=True, type=click.Path(exists=False, file_okay=True, dir_okay=False))
@click.option("--dry-run",  is_flag=True)
@click.option("--verbose", '-v', is_flag=True)
def main(config, config_prefix, population, gid, ref_axis, input_file, template_name, output_file, dry_run, verbose):
    
    utils.config_logging(verbose)
    logger = utils.get_script_logger(os.path.basename(__file__))

    h.load_file("nrngui.hoc")
    h.load_file("import3d.hoc")

    env = Env(config_file=config, config_prefix=config_prefix)
    swc_type_defs = env.SWC_Types

    if not os.path.isfile(output_file):
        io_utils.make_h5types(env, output_file)

    (forest_pop_ranges, _)  = read_population_ranges(output_file)
    (forest_population_start, forest_population_count) = forest_pop_ranges[population]
    forest_population_end = forest_population_start + forest_population_count
    h.load_file(input_file)
    cell = getattr(h, template_name)(0, 0)
    if verbose:
        h.topology()
    tree_dict = export_swc_dict(cell, ref_axis=ref_axis)

    if (gid < forest_population_start) or (gid > forest_population_end):
        gid = forest_population_start
    trees_dict = { gid : tree_dict }

    logger.info(pprint.pformat(trees_dict))

    if not dry_run:
        append_cell_trees(output_file, population, trees_dict)

if __name__ == '__main__':
    main(args=sys.argv[(utils.list_find(lambda x: os.path.basename(x) == os.path.basename(__file__), sys.argv)+1):])


