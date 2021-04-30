import collections, os, sys, traceback, copy, datetime, math, pprint
import numpy as np
#from ca1.neuron_utils import h, d_lambda, default_hoc_sec_lists, default_ordered_sec_types, freq, make_rec, \
#    load_cell_template, HocCellInterface, IzhiCellAttrs, default_izhi_cell_attrs_dict
from ca1.neuron_utils import load_cell_template
from ca1.utils import get_module_logger, map, range, zip, zip_longest, viewitems, read_from_yaml, write_to_yaml, Promise
from neuroh5.io import read_cell_attribute_selection, read_graph_selection, read_tree_selection




def make_neurotree_cell(template_class, gid=0, dataset_path="", neurotree_dict={}):
    """
    :param template_class:
    :param local_id:
    :param gid:
    :param dataset_path:
    :param neurotree_dict:
    :return: hoc cell object
    """
    vx = neurotree_dict['x']
    vy = neurotree_dict['y']
    vz = neurotree_dict['z']
    vradius = neurotree_dict['radius']
    vlayer = neurotree_dict['layer']
    vsection = neurotree_dict['section']
    secnodes = neurotree_dict['section_topology']['nodes']
    vsrc = neurotree_dict['section_topology']['src']
    vdst = neurotree_dict['section_topology']['dst']
    vloc = neurotree_dict['section_topology']['loc']
    swc_type = neurotree_dict['swc_type']
    cell = template_class(gid, dataset_path, secnodes, vlayer, vsrc, vdst, vloc, vx, vy, vz, vradius, swc_type)
    return cell

    
def make_section_graph(neurotree_dict):
    """
    Creates a graph of sections that follows the topological organization of the given neuron.
    :param neurotree_dict:
    :return: NetworkX.DiGraph
    """
    import networkx as nx

    if 'section_topology' in neurotree_dict:
        sec_src = neurotree_dict['section_topology']['src']
        sec_dst = neurotree_dict['section_topology']['dst']
        sec_loc = neurotree_dict['section_topology']['loc']
    else:
        sec_src = neurotree_dict['src']
        sec_dst = neurotree_dict['dst']
        sec_loc = []
        sec_nodes = {}
        pt_sections = neurotree_dict['sections']
        pt_parents = neurotree_dict['parent']
        sec_nodes = make_section_node_dict(neurotree_dict)
        for src, dst in zip_longest(sec_src, sec_dst):
            src_pts = sec_nodes[src]
            dst_pts = sec_nodes[dst]
            dst_parent = pt_parents[dst_pts[0]]
            loc = np.argwhere(src_pts == dst_parent)[0]
            sec_loc.append(loc)
                
    sec_graph = nx.DiGraph()
    for i, j, loc in zip(sec_src, sec_dst, sec_loc):
        sec_graph.add_edge(i, j, loc=loc)

    return sec_graph
