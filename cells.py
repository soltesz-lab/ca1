import collections, os, sys, traceback, copy, datetime, math, pprint
import networkx as nx
import numpy as np
from ca1.neuron_utils import load_cell_template, h, d_lambda, init_nseg, reinit_diam, default_hoc_sec_lists, default_ordered_sec_types, make_rec
from ca1.utils import get_module_logger, map, range, zip, zip_longest, viewitems, read_from_yaml, write_to_yaml, Promise
from neuroh5.io import read_cell_attribute_selection, read_graph_selection, read_tree_selection


# This logger will inherit its settings from the root logger, created in ca1.env
logger = get_module_logger(__name__)

class SectionNode(object):
    __slots__ = []
    
    def __init__(self, section_type, index, section, diam_bounds=None):
        self.section = section
        self.index = index
        self.section_type = section_type
        self.diam_bounds = None
        

    @property
    def sec(self):
        return self.section
    

def make_neurotree_hoc_cell(template_class, gid=0, neurotree_dict={}):
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
    cell = template_class(gid, secnodes, vlayer, vsrc, vdst, vloc, vx, vy, vz, vradius, swc_type)
    
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


class BiophysCell(object):
    """
    A Python wrapper for neuronal cell objects specified in the NEURON language hoc.
    """
    def __init__(self, gid, population_name, hoc_cell=None, neurotree_dict=None, mech_file_path=None, mech_dict=None, env=None):
        """

        :param gid: int
        :param population_name: str
        :param hoc_cell: :class:'h.hocObject': instance of a NEURON cell template
        :param mech_file_path: str (path)
        :param env: :class:'Env'
        """
        
        self.gid = gid
        self.population_name = pop_name
        self.tree = nx.DiGraph()
        self.template_class = None
        
        if env is not None:
            self.template_class = env.template_dict[population_name]
            for sec_type in env.SWC_Types:
                if sec_type not in default_ordered_sec_types:
                    raise AttributeError('Unexpected SWC Type definitions found in Env')
                
        self.sections = {key: [] for key in default_ordered_sec_types}
        self.mech_file_path = mech_file_path
        self.init_mech_dict = dict(mech_dict) if mech_dict is not None else None
        self.mech_dict = dict(mech_dict) if mech_dict is not None else None
        self.spike_detector = None
        self.spike_onset_delay = 0.
        self.hoc_cell = hoc_cell
        if hoc_cell is not None:
            import_morphology_from_hoc(self, hoc_cell)
        elif neurotree_dict is not None:
            hoc_cell = make_neurotree_hoc_cell(template_class, gid, neurotree_dict)
            import_morphology_from_hoc(self, hoc_cell)
        if (mech_dict is None) and (mech_file_path is not None):
            import_mech_dict_from_file(self, self.mech_file_path)
        elif mech_dict is None:
            # Allows for a cell to be created and for a new mech_dict to be constructed programmatically from scratch
            self.init_mech_dict = dict()
            self.mech_dict = dict()
        self.root = None
        sorted_nodes = nx.topological_sort(self.tree)
        if len(sorted_nodes) > 0:
            self.root = sorted_nodes[0]
            
        init_cable(self)
        init_spike_detector(self)

    @property
    def gid(self):
        return self.gid

    @property
    def population_name(self):
        return self.population_name

    @property
    def soma(self):
        return self.sections.soma

    @property
    def axon(self):
        return self.sections.axon

    @property
    def basal(self):
        return self.sections.basal

    @property
    def apical(self):
        return self.sections.apical

    @property
    def trunk(self):
        return self.sections.trunk

    @property
    def tuft(self):
        return self.sections.tuft

    @property
    def spine(self):
        return self.sections.spine

    @property
    def ais(self):
        return self.sections.ais

    @property
    def hillock(self):
        return self.sections.hillock


def get_distance_to_node(cell, node, root=None, loc=None):
    """
    Returns the distance from the given location on the given node to its connection with a root node.
    :param node: int
    :param loc: float
    :return: int or float
    """
    if root is None:
        root = cell.root
        
    length = 0.
    if (node is root) or (root is None) or (node is None):
        return length
    if loc is not None:
        length += loc * node.section.L
    rpath = nx.shortest_path(cell.tree, source=root, target=node).reverse()
    while not rpath.empty():
        node = rpath.pop()
        if not rpath.empty():
            parent = rpath.top()
            e = G.get_edge_data(parent, node)
            loc = e['parent_loc']
            length += loc * parent.section.L
    return length


def insert_section_node(cell, section_type, index, sec):
    node = SectionNode(section_type, index, sec)
    if cell.tree.has_node(node):
        raise RuntimeError(f'insert_section: section index {index} already exists in cell {self.gid}')
    cell.tree.add_node(node)
    cell.sections[section_type].append(node)
    return node
    
def insert_section_tree(cell, sec_list, sec_dict, parent=None, connect_hoc_sections=False):
    sec_stack = []
    for sec in sec_list:
        sec_stack.append((parent, sec))
    while not sec_stack.empty():
        sec_parent, sec = sec_stack.pop()
        sec_info = sec_dict[sec]
        sec_children = sec.children()
        for child in sec_children:
            sec_stack.append((sec, child))
        sec_node = insert_section_node(cell, sec_dict[sec].section_type, sec_dict[sec].index, sec)
        if sec_parent is not None:
            sec_parent_node = sec_dict[parent]
            cell.tree = connect_nodes(cell.tree, sec_parent_node, sec_node,
                                      connect_hoc_sections=connect_hoc_sections)
    

                          
def connect_nodes(tree, parent, child, parent_loc=1., child_loc=0., connect_hoc_sections=False):
    """
    Connects the given section node to a parent node, and if specified, establishes a connection between their associated
    hoc sections.
    :param parent: SectionNode
    :param child: SectionNode
    :param parent_loc: float in [0,1] : connect to this end of the parent hoc section
    :param child_loc: float in [0,1] : connect this end of the child hoc section
    :param connect_hoc_sections: bool
    """
    tree.add_edge(parent, child, parent_loc=parent_loc, child_loc=child_loc)
    if connect_hoc_sections:
        child.section.connect(parent.section, parent_loc, child_loc)
    return tree


def import_morphology_from_hoc(cell, hoc_cell):
    """
    Append sections from an existing instance of a NEURON cell template to a Python cell wrapper.
    :param cell: :class:'BiophysCell'
    :param hoc_cell: :class:'h.hocObject': instance of a NEURON cell template
    """
    sec_info_dict = {}
    root_sec = None
    for sec_type, sec_index_list in viewitems(default_hoc_sec_lists):
        if hasattr(hoc_cell, sec_type) and (getattr(hoc_cell, sec_type) is not None):
            sec_list = list(getattr(hoc_cell, sec_type))
            if hasattr(hoc_cell, sec_index_list):
                sec_indexes = list(getattr(hoc_cell, sec_index_list))
            else:
                raise AttributeError('import_morphology_from_hoc: %s is not an attribute of the hoc cell' %
                                     sec_index_list)
            if sec_type == 'soma':
                root_sec = sec_list[0]
            for sec, index in zip(sec_list, sec_indexes):
                sec_info_dict[sec] = { 'type': sec_type, 'index': int(index) }
    if root_sec:
        insert_section_tree(cell, [root_sec], sec_info_dict)
    else
        raise RuntimeError(f'import_morphology_from_hoc: unable to locate root section')


def import_mech_dict_from_file(cell, mech_file_path=None):
    """
    Imports from a .yaml file a dictionary specifying parameters of NEURON cable properties, density mechanisms, and
    point processes for each type of section in a BiophysCell.
    :param cell: :class:'BiophysCell'
    :param mech_file_path: str (path)
    """
    if mech_file_path is None:
        if cell.mech_file_path is None:
            raise ValueError('import_mech_dict_from_file: missing mech_file_path')
        elif not os.path.isfile(cell.mech_file_path):
            raise IOError('import_mech_dict_from_file: invalid mech_file_path: %s' % cell.mech_file_path)
    elif not os.path.isfile(mech_file_path):
        raise IOError('import_mech_dict_from_file: invalid mech_file_path: %s' % mech_file_path)
    else:
        cell.mech_file_path = mech_file_path
    cell.init_mech_dict = read_from_yaml(cell.mech_file_path)
    cell.mech_dict = copy.deepcopy(cell.init_mech_dict)
    
    

def init_cable(cell, verbose=False):
    for sec_type in cell.nodes:
        for node in cell.nodes[sec_type]:
            reset_cable_by_node(cell, node, verbose=verbose)

            
def reset_cable_by_node(cell, node, verbose=True):
    """
    Consults a dictionary specifying parameters of NEURON cable properties such as axial resistance ('Ra'),
    membrane specific capacitance ('cm'), and a spatial resolution parameter to specify the number of separate
    segments per section in a BiophysCell
    :param cell: :class:'BiophysCell'
    :param node_index: int
    :param verbose: bool
    """
    sec_type = node.section_type
    if sec_type in cell.mech_dict and 'cable' in cell.mech_dict[sec_type]:
        mech_content = cell.mech_dict[sec_type]['cable']
        if mech_content is not None:
            update_mechanisms_by_node(cell, node, 'cable', mech_content, verbose=verbose)
    else:
        init_nseg(node.section, verbose=verbose)
        if hasattr(node, 'diam_bounds'):
            reinit_diam(node.section, getattr(node, 'diam_bounds'))

        
def connect2target(cell, sec, loc=1., param='_ref_v', delay=None, weight=None, threshold=None, target=None):
    """
    Converts analog voltage in the specified section to digital spike output. Initializes and returns an h.NetCon
    object with voltage as a reference parameter connected to the specified target.
    :param cell: :class:'BiophysCell'
    :param sec: :class:'h.Section'
    :param loc: float
    :param param: str
    :param delay: float
    :param weight: float
    :param threshold: float
    :param target: object that can receive spikes
    :return: :class:'h.NetCon'
    """
    if cell.spike_detector is not None:
        if delay is None:
            delay = cell.spike_detector.delay
        if weight is None:
            weight = cell.spike_detector.weight[0]
        if threshold is None:
            threshold = cell.spike_detector.threshold
    else:
        if delay is None:
            delay = 0.
        if weight is None:
            weight = 1.
        if threshold is None:
            threshold = -30.
    ps = getattr(sec(loc), param)
    this_netcon = h.NetCon(ps, target, sec=sec)
    this_netcon.delay = delay
    this_netcon.weight[0] = weight
    this_netcon.threshold = threshold
    return this_netcon
            

def init_spike_detector(cell, node=None, distance=100., threshold=-30, delay=0.05, onset_delay=0., loc=0.5):
    """
    Initializes the spike detector in the given cell according to the
    given arguments or a spike detector configuration of the mechanism
    dictionary of the cell, if one exists.

    :param cell: :class:'BiophysCell'
    :param node: :class:'SectionNode'
    :param distance: float
    :param threshold: float
    :param delay: float
    :param onset_delay: float
    :param loc: float
    """
    if cell.mech_dict is not None:
        if 'spike detector' in cell.mech_dict:
            config = cell.mech_dict['spike detector']
            node = getattr(cell, config['section'])[0]
            loc = config['loc']
            distance = config['distance']
            threshold = config['threshold']
            delay = config['delay']
            onset_delay = config['onset delay']

    if node is None:
        if cell.axon:
            for node in cell.axon:
                sec_seg_locs = [seg.x for seg in node.sec]
                for loc in sec_seg_locs:
                    if get_distance_to_node(cell, cell.root, node, loc=loc) >= distance:
                        break
                else:
                    continue
                break
            else:
                node = cell.axon[-1]
                loc = 1.
        elif cell.ais:
            node = cell.ais[0]
        elif cell.soma:
            node = cell.soma[0]
        else:
            raise RuntimeError('init_spike_detector: cell has neither soma nor axon compartment')

    cell.spike_detector = connect2target(cell, node.section, loc=loc, delay=delay, threshold=threshold)

    cell.onset_delay = onset_delay
            
    return cell.spike_detector


def update_mechanism_by_node(cell, node, mech_name, mech_content, verbose=True):
    """
    This method loops through all the parameters for a single mechanism specified in the mechanism dictionary and
    calls apply_mech_rules to interpret the rules and set the values for the given node.
    :param cell: :class:'BiophysCell'
    :param node: :class:'SectionNode'
    :param mech_name: str
    :param mech_content: dict
    :param verbose: bool
    """
    if mech_content is not None:
        for param_name in mech_content:
            # accommodate either a dict, or a list of dicts specifying multiple location constraints for
            # a single parameter
            if isinstance(mech_content[param_name], dict):
                apply_mech_rules(cell, node, mech_name, param_name, mech_content[param_name], verbose=verbose)
            elif isinstance(mech_content[param_name], list):
                for mech_content_entry in mech_content[param_name]:
                    apply_mech_rules(cell, node, mech_name, param_name, mech_content_entry, verbose=verbose)
    else:
        node.section.insert(mech_name)
