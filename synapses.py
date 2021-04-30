import sys, collections, copy, itertools, math, pprint, time, traceback
from functools import reduce
from collections import defaultdict
import numpy as np
import scipy
import scipy.optimize as opt
from neuroh5.io import write_cell_attributes
from ca1.cells import make_section_graph
#from ca1.cells import get_distance_to_node, get_donor, get_mech_rules_dict, get_param_val_by_distance, \
#    import_mech_dict_from_file, make_section_graph, custom_filter_if_terminal, \
#    custom_filter_modify_slope_if_terminal, custom_filter_by_branch_order
#from ca1.neuron_utils import h, default_ordered_sec_types, mknetcon, mknetcon_vecstim
from ca1.utils import ExprClosure, Promise, NamedTupleWithDocstring, get_module_logger, generator_ifempty, map, range, str, \
     viewitems, viewkeys, zip, zip_longest, partitionn, rejection_sampling

# This logger will inherit its settings from the root logger, created in dentate.env
logger = get_module_logger(__name__)

def get_node_attribute(name, content, sec, secnodes, x=None):
    """

    :param name:
    :param content:
    :param sec:
    :param secnodes:
    :param x:
    :return:
    """
    if name in content:
        if x is None:
            return content[name]
        elif sec.n3d() == 0:
            return content[name][0]
        else:
            prev = None
            for i in range(sec.n3d()):
                pos = (sec.arc3d(i) / sec.L)
                if pos >= x:
                    if (prev is None) or (abs(pos - x) < abs(prev - x)):
                        return content[name][secnodes[i]]
                    else:
                        return content[name][secnodes[i - 1]]
                else:
                    prev = pos
    else:
        return None


def synapse_seg_density(syn_type_dict, layer_dict, layer_density_dicts, seg_dict, ran, neurotree_dict=None):
    """
    Computes per-segment density of synapse placement.
    :param syn_type_dict:
    :param layer_dict:
    :param layer_density_dicts:
    :param seg_dict:
    :param ran:
    :param neurotree_dict:
    :return:
    """
    segdensity_dict = {}
    layers_dict = {}

    if neurotree_dict is not None:
        secnodes_dict = neurotree_dict['section_topology']['nodes']
    else:
        secnodes_dict = None
    for (syn_type_label, layer_density_dict) in viewitems(layer_density_dicts):
        syn_type = syn_type_dict[syn_type_label]
        rans = {}
        for (layer_label, density_dict) in viewitems(layer_density_dict):
            if layer_label == 'default':
                layer = layer_label
            else:
                layer = int(layer_dict[layer_label])
            rans[layer] = ran
        segdensity = defaultdict(list)
        layers = defaultdict(list)
        total_seg_density = 0.
        for sec_index, seg_list in viewitems(seg_dict):
            for seg in seg_list:
                L = seg.sec.L
                nseg = seg.sec.nseg
                if neurotree_dict is not None:
                    secnodes = secnodes_dict[sec_index]
                    layer = get_node_attribute('layer', neurotree_dict, seg.sec, secnodes, seg.x)
                else:
                    layer = -1
                layers[sec_index].append(layer)

                this_ran = None

                if layer > -1:
                    if layer in rans:
                        this_ran = rans[layer]
                    elif 'default' in rans:
                        this_ran = rans['default']
                    else:
                        this_ran = None
                elif 'default' in rans:
                    this_ran = rans['default']
                else:
                    this_ran = None
                if this_ran is not None:
                    while True:
                        dens = this_ran.normal(density_dict['mean'], density_dict['variance'])
                        if dens > 0.0:
                            break
                else:
                    dens = 0.
                total_seg_density += dens
                segdensity[sec_index].append(dens)

        if total_seg_density < 1e-6:
            logger.warning("sections with zero %s synapse density: %s; rans: %s; density_dict: %s; morphology: %s" % (
            syn_type_label, str(segdensity), str(rans), str(density_dict), str(neurotree_dict)))

        segdensity_dict[syn_type] = segdensity

        layers_dict[syn_type] = layers
    return (segdensity_dict, layers_dict)


def synapse_seg_counts(syn_type_dict, layer_dict, layer_density_dicts, sec_index_dict, seg_dict, ran,
                       neurotree_dict=None):
    """
    Computes per-segment relative counts of synapse placement.
    :param syn_type_dict:
    :param layer_dict:
    :param layer_density_dicts:
    :param sec_index_dict:
    :param seg_dict:
    :param seed:
    :param neurotree_dict:
    :return:
    """
    segcounts_dict = {}
    layers_dict = {}
    segcount_total = 0
    if neurotree_dict is not None:
        secnodes_dict = neurotree_dict['section_topology']['nodes']
    else:
        secnodes_dict = None
    for (syn_type_label, layer_density_dict) in viewitems(layer_density_dicts):
        syn_type = syn_type_dict[syn_type_label]
        rans = {}
        for (layer_label, density_dict) in viewitems(layer_density_dict):
            if layer_label == 'default':
                layer = layer_label
            else:
                layer = layer_dict[layer_label]

            rans[layer] = ran
        segcounts = []
        layers = []
        for sec_index, seg_list in viewitems(seg_dict):
            for seg in seg_list:
                L = seg.sec.L
                nseg = seg.sec.nseg
                if neurotree_dict is not None:
                    secnodes = secnodes_dict[sec_index]
                    layer = get_node_attribute('layer', neurotree_dict, seg.sec, secnodes, seg.x)
                else:
                    layer = -1
                layers.append(layer)

                ran = None

                if layer > -1:
                    if layer in rans:
                        ran = rans[layer]
                    elif 'default' in rans:
                        ran = rans['default']
                    else:
                        ran = None
                elif 'default' in rans:
                    ran = rans['default']
                else:
                    ran = None
                if ran is not None:
                    l = (L / nseg)
                    dens = ran.normal(density_dict['mean'], density_dict['variance'])
                    rc = dens * l
                    segcount_total += rc
                    segcounts.append(rc)
                else:
                    segcounts.append(0)

            segcounts_dict[syn_type] = segcounts
            layers_dict[syn_type] = layers
    return (segcounts_dict, segcount_total, layers_dict)



def distribute_uniform_synapses(density_seed, syn_type_dict, swc_type_dict, layer_dict, sec_layer_density_dict,
                                neurotree_dict, cell_sec_dict, cell_secidx_dict):
    """
    Computes uniformly-spaced synapse locations.
    :param density_seed:
    :param syn_type_dict:
    :param swc_type_dict:
    :param layer_dict:
    :param sec_layer_density_dict:
    :param neurotree_dict:
    :param sec_dict:
    :param secidx_dict:
    :return:
    """
    syn_ids = []
    syn_locs = []
    syn_secs = []
    syn_layers = []
    syn_types = []
    swc_types = []
    syn_index = 0

    r = np.random.RandomState()
    local_random.seed(int(seed))

    segcounts_per_sec = {}
    for (sec_name, layer_density_dict) in viewitems(sec_layer_density_dict):
        sec_index_dict = cell_secidx_dict[sec_name]
        swc_type = swc_type_dict[sec_name]
        seg_list = []
        L_total = 0
        (seclst, maxdist) = cell_sec_dict[sec_name]
        secidxlst = cell_secidx_dict[sec_name]
        sec_dict = {int(idx): sec for sec, idx in zip(seclst, secidxlst)}
        seg_dict = {}
        for (sec_index, sec) in viewitems(sec_dict):
            seg_list = []
            if maxdist is None:
                for seg in sec:
                    if seg.x < 1.0 and seg.x > 0.0:
                        seg_list.append(seg)
            else:
                for seg in sec:
                    if seg.x < 1.0 and seg.x > 0.0 and ((L_total + sec.L * seg.x) <= maxdist):
                        seg_list.append(seg)
            L_total += sec.L
            seg_dict[sec_index] = seg_list
        segcounts_dict, total, layers_dict = \
            synapse_seg_counts(syn_type_dict, layer_dict, layer_density_dict, \
                               sec_index_dict=sec_index_dict, seg_dict=seg_dict, ran=r, \
                               neurotree_dict=neurotree_dict)
        segcounts_per_sec[sec_name] = segcounts_dict
        sample_size = total
        for (syn_type_label, _) in viewitems(layer_density_dict):
            syn_type = syn_type_dict[syn_type_label]
            segcounts = segcounts_dict[syn_type]
            layers = layers_dict[syn_type]
            for sec_index, seg_list in viewitems(seg_dict):
                for seg, layer, seg_count in zip(seg_list, layers, segcounts):
                    seg_start = seg.x - (0.5 / seg.sec.nseg)
                    seg_end = seg.x + (0.5 / seg.sec.nseg)
                    seg_range = seg_end - seg_start
                    int_seg_count = math.floor(seg_count)
                    syn_count = 0
                    while syn_count < int_seg_count:
                        syn_loc = seg_start + seg_range * (syn_count + 1) / math.ceil(seg_count)
                        assert ((syn_loc <= 1) & (syn_loc >= 0))
                        if syn_loc < 1.0:
                            syn_locs.append(syn_loc)
                            syn_ids.append(syn_index)
                            syn_secs.append(sec_index_dict[seg.sec])
                            syn_layers.append(layer)
                            syn_types.append(syn_type)
                            swc_types.append(swc_type)
                            syn_index += 1
                            syn_count += 1

    assert (len(syn_ids) > 0)
    syn_dict = {'syn_ids': np.asarray(syn_ids, dtype='uint32'),
                'syn_locs': np.asarray(syn_locs, dtype='float32'),
                'syn_secs': np.asarray(syn_secs, dtype='uint32'),
                'syn_layers': np.asarray(syn_layers, dtype='int8'),
                'syn_types': np.asarray(syn_types, dtype='uint8'),
                'swc_types': np.asarray(swc_types, dtype='uint8')}

    return (syn_dict, segcounts_per_sec)


def distribute_poisson_synapses(density_seed, syn_type_dict, swc_type_dict, layer_dict, sec_layer_density_dict,
                                neurotree_dict, cell_sec_dict, cell_secidx_dict):
    """
    Computes synapse locations distributed according to a Poisson distribution.
    :param density_seed:
    :param syn_type_dict:
    :param swc_type_dict:
    :param layer_dict:
    :param sec_layer_density_dict:
    :param neurotree_dict:
    :param cell_sec_dict:
    :param cell_secidx_dict:
    :param verbose:
    :return:
    """
    import networkx as nx

    syn_ids = []
    syn_locs = []
    syn_secs = []
    syn_layers = []
    syn_types = []
    swc_types = []
    syn_index = 0

    sec_graph = make_section_graph(neurotree_dict)

    debug_flag = False
    secnodes_dict = neurotree_dict['section_topology']['nodes']
    for sec, secnodes in viewitems(secnodes_dict):
        if len(secnodes) < 2:
            debug_flag = True

    if debug_flag:
        logger.debug('sec_graph: %s' % str(list(sec_graph.edges)))
        logger.debug('neurotree_dict: %s' % str(neurotree_dict))

    seg_density_per_sec = {}
    r = np.random.RandomState()
    r.seed(int(density_seed))
    for (sec_name, layer_density_dict) in viewitems(sec_layer_density_dict):

        swc_type = swc_type_dict[sec_name]
        seg_dict = {}
        L_total = 0

        (seclst, maxdist) = cell_sec_dict[sec_name]
        secidxlst = cell_secidx_dict[sec_name]
        sec_dict = {int(idx): sec for sec, idx in zip(seclst, secidxlst)}
        if len(sec_dict) > 1:
            sec_subgraph = sec_graph.subgraph(list(sec_dict.keys()))
            if len(sec_subgraph.edges()) > 0:
                sec_roots = [n for n, d in sec_subgraph.in_degree() if d == 0]
                sec_edges = []
                for sec_root in sec_roots:
                    sec_edges.append(list(nx.dfs_edges(sec_subgraph, sec_root)))
                    sec_edges.append([(None, sec_root)])
                sec_edges = [val for sublist in sec_edges for val in sublist]
            else:
                sec_edges = [(None, idx) for idx in list(sec_dict.keys())]
        else:
            sec_edges = [(None, idx) for idx in list(sec_dict.keys())]
        for sec_index, sec in viewitems(sec_dict):
            seg_list = []
            if maxdist is None:
                for seg in sec:
                    if seg.x < 1.0 and seg.x > 0.0:
                        seg_list.append(seg)
            else:
                for seg in sec:
                    if seg.x < 1.0 and seg.x > 0.0 and ((L_total + sec.L * seg.x) <= maxdist):
                        seg_list.append(seg)
            seg_dict[sec_index] = seg_list
            L_total += sec.L
        seg_density_dict, layers_dict = \
            synapse_seg_density(syn_type_dict, layer_dict, \
                                layer_density_dict, \
                                seg_dict, r, \
                                neurotree_dict=neurotree_dict)
        seg_density_per_sec[sec_name] = seg_density_dict
        for (syn_type_label, _) in viewitems(layer_density_dict):
            syn_type = syn_type_dict[syn_type_label]
            seg_density = seg_density_dict[syn_type]
            layers = layers_dict[syn_type]
            end_distance = {}
            for sec_parent, sec_index in sec_edges:
                seg_list = seg_dict[sec_index]
                sec_seg_layers = layers[sec_index]
                sec_seg_density = seg_density[sec_index]
                start_seg = seg_list[0]
                interval = 0.
                syn_loc = 0.
                for seg, layer, density in zip(seg_list, sec_seg_layers, sec_seg_density):
                    seg_start = seg.x - (0.5 / seg.sec.nseg)
                    seg_end = seg.x + (0.5 / seg.sec.nseg)
                    L = seg.sec.L
                    L_seg_start = seg_start * L
                    L_seg_end = seg_end * L
                    if density > 0.:
                        beta = 1. / density
                        if interval > 0.:
                            sample = r.exponential(beta)
                        else:
                            while True:
                                sample = r.exponential(beta)
                                if (sample >= L_seg_start) and (sample < L_seg_end):
                                    break
                        interval += sample
                        while interval < L_seg_end:
                            if interval >= L_seg_start:
                                syn_loc = (interval / L)
                                assert ((syn_loc <= 1) and (syn_loc >= seg_start))
                                if syn_loc < 1.0:
                                    syn_locs.append(syn_loc)
                                    syn_ids.append(syn_index)
                                    syn_secs.append(sec_index)
                                    syn_layers.append(layer)
                                    syn_types.append(syn_type)
                                    swc_types.append(swc_type)
                                    syn_index += 1
                            interval += r.exponential(beta)
                    else:
                        interval = seg_end * L
                end_distance[sec_index] = (1.0 - syn_loc) * L

    assert (len(syn_ids) > 0)
    syn_dict = {'syn_ids': np.asarray(syn_ids, dtype='uint32'),
                'syn_locs': np.asarray(syn_locs, dtype='float32'),
                'syn_secs': np.asarray(syn_secs, dtype='uint32'),
                'syn_layers': np.asarray(syn_layers, dtype='int8'),
                'syn_types': np.asarray(syn_types, dtype='uint8'),
                'swc_types': np.asarray(swc_types, dtype='uint8')}

    return (syn_dict, seg_density_per_sec)

