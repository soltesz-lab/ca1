import numbers, os, copy, pprint, sys
from collections import defaultdict
from scipy import interpolate, signal
import numpy as np
from mpi4py import MPI
from neuroh5.io import NeuroH5ProjectionGen, bcast_cell_attributes, read_cell_attributes, read_population_names, \
    read_population_ranges, read_projection_names, read_tree_selection
import h5py
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.animation import FuncAnimation, writers
from matplotlib.colors import BoundaryNorm
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import ca1
from ca1.env import Env
from ca1.utils import get_module_logger, Struct, viewitems, make_geometric_graph, zip_longest, signal_psd
from ca1.io_utils import get_h5py_attr, set_h5py_attr
from ca1 import spikedata

try:
    from neural_geometry.geometry import measure_distance_extents, get_total_extents
except ImportError as e:
    print(('neural_geometry.plot: problem importing module required by dentate.geometry:', e))

# This logger will inherit its settings from the root logger, created in ca1.env
logger = get_module_logger(__name__)

# Default figure configuration
default_fig_options = Struct(figFormat='png', lw=2, figSize=(15,8), fontSize=14, saveFig=None, showFig=True,
                             colormap='jet', saveFigDir=None)

dflt_colors = ["#009BFF", "#E85EBE", "#00FF00", "#0000FF", "#FF0000", "#01FFFE", "#FFA6FE", 
              "#FFDB66", "#006401", "#010067", "#95003A", "#007DB5", "#FF00F6", "#FFEEE8", "#774D00",
              "#90FB92", "#0076FF", "#D5FF00", "#FF937E", "#6A826C", "#FF029D", "#FE8900", "#7A4782",
              "#7E2DD2", "#85A900", "#FF0056", "#A42400", "#00AE7E", "#683D3B", "#BDC6FF", "#263400",
              "#BDD393", "#00B917", "#9E008E", "#001544", "#C28C9F", "#FF74A3", "#01D0FF", "#004754",
              "#E56FFE", "#788231", "#0E4CA1", "#91D0CB", "#BE9970", "#968AE8", "#BB8800", "#43002C",
              "#DEFF74", "#00FFC6", "#FFE502", "#620E00", "#008F9C", "#98FF52", "#7544B1", "#B500FF",
              "#00FF78", "#FF6E41", "#005F39", "#6B6882", "#5FAD4E", "#A75740", "#A5FFD2", "#FFB167"]

rainbow_colors = ["#9400D3", "#4B0082", "#00FF00", "#FFFF00", "#FF7F00", "#FF0000"]

raster_colors = ['#8dd3c7', '#ffed6f', '#bebada', '#fb8072', '#80b1d3', '#fdb462',
                    '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5']


def hex2rgb(hexcode):
    if hasattr(hexcode, 'decode'):
        return tuple([ float(b)/255.0 for b in map(ord,hexcode[1:].decode('hex')) ])
    else:
        import codecs
        bhexcode = bytes(hexcode[1:], 'utf-8')
        return tuple([ float(b)/255.0 for b in codecs.decode(bhexcode, 'hex') ])

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.size'] = 14.
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['text.usetex'] = False


def show_figure():
    try:
        plt.show(block=False)
    except:
        plt.show()

def close_figure(fig):
    plt.close(fig)


def save_figure(file_name_prefix, fig=None, **kwargs):
    """

    :param file_name_prefix:
    :param fig: :class:'plt.Figure'
    :param kwargs: dict
    """
    fig_options = copy.copy(default_fig_options)
    fig_options.update(kwargs)
    fig_file_path = '%s.%s' % (file_name_prefix, fig_options.figFormat)
    if fig_options.saveFigDir is not None:
        fig_file_path = '%s/%s' % (fig_options.saveFigDir, fig_file_path)
    if fig is not None:
        fig.savefig(fig_file_path)
    else:
        plt.savefig(fig_file_path)


def plot_graph(x, y, z, start_idx, end_idx, edge_scalars=None, edge_color=None, **kwargs):
    """ 
    Shows graph edges using Mayavi

    Parameters
    -----------
        x: ndarray
            x coordinates of the points
        y: ndarray
            y coordinates of the points
        z: ndarray
            z coordinates of the points
        edge_scalars: ndarray, optional
            optional data to give the color of the edges.
        kwargs:
            extra keyword arguments are passed to quiver3d.
    """
    from mayavi import mlab
    if edge_color is not None:
        kwargs['color'] = edge_color
    vec = mlab.quiver3d(x[start_idx],
                        y[start_idx],
                        z[start_idx],
                        x[end_idx] - x[start_idx],
                        y[end_idx] - y[start_idx],
                        z[end_idx] - z[start_idx],
                        scalars=edge_scalars,
                        scale_factor=1,
                        mode='2ddash',
                        **kwargs)
    b = mlab.points3d(x[0],y[0],z[0],
                      mode='cone',
                      scale_factor=10,
                      **kwargs)
    if edge_scalars is not None:
        vec.glyph.color_mode = 'color_by_scalar'
        cb = mlab.colorbar(vec, label_fmt='%.1f')
        cb.label_text_property.font_size=14
    return vec


def plot_spatial_bin_graph(graph_dict, **kwargs):
    
    import hiveplot as hv
    import networkx as nx
    
    edge_dflt_colors = ['red','crimson','coral','purple']
    
    fig_options = copy.copy(default_fig_options)
    fig_options.update(kwargs)

    label = graph_dict['label']
    GU = graph_dict['U graph']

    destination = graph_dict['destination']
    sources = graph_dict['sources']

    nodes = {}
    nodes[destination] = [(s,d) for s, d in GU.nodes() if s == destination]
    for source in sources:
        nodes[source] = [(s,d) for s, d in GU.nodes() if s == source]

    snodes = {}
    for group, nodelist in viewitems(nodes):
        snodes[group] = sorted(nodelist)

    edges = {}
    for source in sources:
        edges[source] = [(u,v,d) for u,v,d in GU.edges(data=True) if v[0] == source]

    nodes_cmap = dict()
    nodes_cmap[destination] = 'blue'
    for i, source in enumerate(sources):
        nodes_cmap[source] = raster_colors[i]

    edges_cmap = dict()
    for i, source in enumerate(sources):
        edges_cmap[source] = dflt_colors[i]

    hvpl = hv.HivePlot(snodes, edges, nodes_cmap, edges_cmap)
    hvpl.draw()

    filename = '%s.%s' % (label, fig_options.figFormat)
    plt.savefig(filename)
    


def plot_coordinates(coords_path, population, namespace, index = 0, graph_type = 'scatter', bin_size = 0.01, xyz = False, **kwargs):
    """
    Plot coordinates

    :param coords_path:
    :param namespace: 
    :param population: 

    """
    fig_options = copy.copy(default_fig_options)
    fig_options.update(kwargs)
        
    soma_coords = read_cell_attributes(coords_path, population, namespace=namespace)
    
        
    fig = plt.figure(1, figsize=plt.figaspect(1.) * 2.)
    ax = plt.gca()

    coord_U = {}
    coord_V = {}
    if xyz:
        for k,v in soma_coords:
            coord_U[k] = v['X Coordinate'][index]
            coord_V[k] = v['Y Coordinate'][index]
    else:
        for k,v in soma_coords:
            coord_U[k] = v['U Coordinate'][index]
            coord_V[k] = v['V Coordinate'][index]
    
    coord_U_array = np.asarray([coord_U[k] for k in sorted(coord_U.keys())])
    coord_V_array = np.asarray([coord_V[k] for k in sorted(coord_V.keys())])

    x_min = np.min(coord_U_array)
    x_max = np.max(coord_U_array)
    y_min = np.min(coord_V_array)
    y_max = np.max(coord_V_array)

    dx = int((x_max - x_min) / bin_size)
    dy = int((y_max - y_min) / bin_size)

    if graph_type == 'scatter':
        ax.scatter(coord_U_array, coord_V_array, alpha=0.1, linewidth=0)
        ax.axis([x_min, x_max, y_min, y_max])
    elif graph_type == 'histogram2d':
        (H, xedges, yedges) = np.histogram2d(coord_U_array, coord_V_array, bins=[dx, dy])
        X, Y = np.meshgrid(xedges, yedges)
        Hint = H[:-1, :-1]
        levels = MaxNLocator(nbins=25).tick_values(Hint.min(), Hint.max())
        cmap = plt.get_cmap('jet')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        p = ax.contourf(X[:-1,:-1] + (bin_size / 2), Y[:-1,:-1]+(bin_size / 2), H.T, levels=levels, cmap=cmap)
        fig.colorbar(p, ax=ax, shrink=0.5, aspect=20)
    else:
        raise ValueError('Unknown graph type %s' % graph_type)

    if xyz:
        ax.set_xlabel('X coordinate (um)', fontsize=fig_options.fontSize)
        ax.set_ylabel('Y coordinate (um)', fontsize=fig_options.fontSize)
    else:
        ax.set_xlabel('U coordinate (septal - temporal)', fontsize=fig_options.fontSize)
        ax.set_ylabel('V coordinate (supra - infrapyramidal)', fontsize=fig_options.fontSize)
        
    ax.set_title('Coordinate distribution for population: %s' % (population), fontsize=fig_options.fontSize)
    
    if fig_options.saveFig:
        if isinstance(fig_options.saveFig, basestring):
            filename = fig_options.saveFig
        else:
            filename = population+' Coordinates.%s' % fig_options.figFormat
            plt.savefig(filename)

    if fig_options.showFig:
        show_figure()
    
    return ax



def plot_coords_in_volume(populations, coords_path, coords_namespace, config, scale=25., subpopulation=-1, subvol=False, verbose=False, mayavi=False):
    
    env = Env(config_file=config)

    rotate = env.geometry['Parametric Surface']['Rotation']
    layer_extents = env.geometry['Parametric Surface']['Layer Extents']
    rotate = env.geometry['Parametric Surface']['Rotation']

    (extent_u, extent_v, extent_l) = get_total_extents(layer_extents)

    logger.info('Reading coordinates...')

    pop_min_extent = None
    pop_max_extent = None

    xcoords = []
    ycoords = []
    zcoords = []
    cmap = cm.get_cmap('Dark2')
    cmap_range = np.linspace(0,1,num=len(populations))

    colors = []
    for (pop_id, population) in enumerate(populations):
        coords = read_cell_attributes(coords_path, population, namespace=coords_namespace)

        count = 0
        cxcoords = []
        cycoords = []
        czcoords = []
        for (k,v) in coords:
            count += 1
            cxcoords.append(v['X Coordinate'][0])
            cycoords.append(v['Y Coordinate'][0])
            czcoords.append(v['Z Coordinate'][0])
        if subpopulation > -1 and count > subpopulation:
            ridxs  = np.random.choice(np.arange(len(cxcoords)), replace=False, size=subpopulation)
            cxcoords = list(np.asarray(cxcoords)[ridxs])
            cycoords = list(np.asarray(cycoords)[ridxs])
            czcoords = list(np.asarray(czcoords)[ridxs])

        colors += [cmap(cmap_range[pop_id]) for _ in range(len(cxcoords))]
        xcoords += cxcoords
        ycoords += cycoords
        zcoords += czcoords
        logger.info(f'Read {count} coordinates...')
        
        pop_distribution = env.geometry['Cell Distribution'][population]
        pop_layers = []
        for layer in pop_distribution:
            num_layer = pop_distribution[layer]
            if num_layer > 0:
                pop_layers.append(layer)
            
                if pop_min_extent is None:
                    pop_min_extent = np.asarray(layer_extents[layer][0])
                else:
                    pop_min_extent = np.minimum(pop_min_extent, np.asarray(layer_extents[layer][0]))

                if pop_max_extent is None:
                    pop_max_extent = np.asarray(layer_extents[layer][1])
                else:
                    pop_max_extent = np.maximum(pop_min_extent, np.asarray(layer_extents[layer][1]))


    pts = np.concatenate((np.asarray(xcoords).reshape(-1,1), \
                          np.asarray(ycoords).reshape(-1,1), \
                          np.asarray(zcoords).reshape(-1,1)),axis=1)

    if mayavi:
        from mayavi import mlab
    else:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

    
    logger.info('Plotting coordinates...')
    if mayavi: 
        mlab.points3d(*pts.T, color=(1, 1, 0), scale_factor=scale)
    else:
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection='3d')
        ax.scatter(*pts.T, c=colors, s=int(scale))
        
    logger.info('Constructing volume...')
    from ca1.CA1_volume import make_CA1_volume

    if subvol:
        subvol = make_CA1_volume ((pop_min_extent[0], pop_max_extent[0]), \
                                (pop_min_extent[1], pop_max_extent[1]), \
                                (pop_min_extent[2], pop_max_extent[2]), \
                                resolution=[3, 3, 3], \
                                rotate=rotate)
    else:
        vol = make_CA1_volume ((extent_u[0], extent_u[1]),
                            (extent_v[0], extent_v[1]),
                            (extent_l[0], extent_l[1]),
                            resolution=[3, 3, 3],
                            rotate=rotate)

    logger.info('Plotting volume...')

    if subvol:
        if mayavi:
            subvol.mplot_surface(color=(0, 0.4, 0), opacity=0.33)
        else:
            subvol.mplot_surface(color='k', alpha=0.33, figax=[fig, ax])
    else:
        if mayavi:
            vol.mplot_surface(color=(0, 1, 0), opacity=0.33)
        else:
            vol.mplot_surface(color='k', alpha=0.33, figax=[fig, ax])
    if mayavi:
        mlab.show()
    else:
        ax.view_init(-90,0)
        plt.show()


        
def plot_cell_tree(population, gid, tree_dict, line_width=1., sample=0.05, color_edge_scalars=True, mst=False, conn_loc=True):
    
    import networkx as nx
    from mayavi import mlab
    
    mlab.figure(bgcolor=(0,0,0))

    xcoords = tree_dict['x']
    ycoords = tree_dict['y']
    zcoords = tree_dict['z']
    swc_type = tree_dict['swc_type']
    layer    = tree_dict['layer']
    secnodes = tree_dict['section_topology']['nodes']
    src      = tree_dict['section_topology']['src']
    dst      = tree_dict['section_topology']['dst']
    loc      = tree_dict['section_topology']['loc']
    
    x = xcoords.reshape(-1,)
    y = ycoords.reshape(-1,)
    z = zcoords.reshape(-1,)

    edges = []
    for sec, nodes in viewitems(secnodes):
        for i in range(1, len(nodes)):
            srcnode = nodes[i-1]
            dstnode = nodes[i]
            edges.append((srcnode, dstnode))

    loc_x = []
    loc_y = []
    loc_z = []
    for (s,d,l) in zip(src,dst,loc):
        srcnode = secnodes[s][l]
        dstnode = secnodes[d][0]
        edges.append((srcnode, dstnode))
        loc_x.append(x[srcnode])
        loc_y.append(y[srcnode])
        loc_z.append(z[srcnode])

    conn_loc_x = np.asarray(loc_x, dtype=np.float64)
    conn_loc_y = np.asarray(loc_y, dtype=np.float64)
    conn_loc_z = np.asarray(loc_z, dtype=np.float64)
        
    # Make a NetworkX graph out of our point and edge data
    g = make_geometric_graph(x, y, z, edges)

    edges = g.edges
    # Compute minimum spanning tree using networkx
    # nx.mst returns an edge generator
    if mst:
        edges = nx.minimum_spanning_tree(g).edges(data=True)

    edge_array = np.array(list(edges)).T
    start_idx = edge_array[0, :]
    end_idx = edge_array[1, :]

    
    start_idx = start_idx.astype(np.int)
    end_idx   = end_idx.astype(np.int)
    if color_edge_scalars:
        edge_scalars = z[start_idx]
        edge_color = None
    else:
        edge_scalars = None
        edge_color = hex2rgb(rainbow_colors[gid%len(rainbow_colors)])
        
    fig = mlab.gcf()
    
    # Plot this with Mayavi
    g = plot_graph(x, y, z, start_idx, end_idx, edge_scalars=edge_scalars, edge_color=edge_color, \
                   opacity=0.8, colormap='summer', line_width=line_width, figure=fig)

    if conn_loc:
       conn_pts = mlab.points3d(conn_loc_x, conn_loc_y, conn_loc_z, figure=fig,
                                mode='2dcross', colormap='copper', scale_factor=10)

    
    fig.scene.x_plus_view()
    mlab.savefig('%s_%d_cell_tree.x3d' % (population, gid), figure=fig, magnification=10)
    mlab.show()
        

    
## Plot spike raster
def plot_spike_raster (input_path, namespace_id, include = ['eachPop'], time_range = None, time_variable='t', max_spikes = int(1e6), labels = 'legend', pop_rates = True, spike_hist = None, spike_hist_bin = 5, include_artificial=True, marker='.', **kwargs):
    ''' 
    Raster plot of network spike times. Returns the figure handle.

    input_path: file with spike data
    namespace_id: attribute namespace for spike events
    time_range ([start:stop]): Time range of spikes shown; if None shows all (default: None)
    time_variable: Name of variable containing spike times (default: 't')
    max_spikes (int): maximum number of spikes that will be plotted  (default: 1e6)
    labels = ('legend', 'overlay'): Show population labels in a legend or overlayed on one side of raster (default: 'legend')
    pop_rates = (True|False): Include population rates (default: False)
    spike_hist (None|'overlay'|'subplot'): overlay line over raster showing spike histogram (spikes/bin) (default: False)
    spike_hist_bin (int): Size of bin in ms to use for histogram (default: 5)
    marker (char): Marker for each spike (default: '|')
    '''

    fig_options = copy.copy(default_fig_options)
    fig_options.update(kwargs)

    (population_ranges, N) = read_population_ranges(input_path)
    population_names  = read_population_names(input_path)

    total_num_cells = 0
    pop_num_cells = {}
    pop_start_inds = {}
    for k in population_names:
        pop_start_inds[k] = population_ranges[k][0]
        pop_num_cells[k] = population_ranges[k][1]
        total_num_cells += population_ranges[k][1]

    include = list(include)
    # Replace 'eachPop' with list of populations
    if 'eachPop' in include: 
        include.remove('eachPop')
        for pop in population_names:
            include.append(pop)
            
    # sort according to start index        
    include.sort(key=lambda x: pop_start_inds[x])
    
    spkdata = spikedata.read_spike_events (input_path, include, namespace_id,
                                           include_artificial=include_artificial,
                                           spike_train_attr_name=time_variable,
                                           time_range=time_range)

    spkpoplst        = spkdata['spkpoplst']
    spkindlst        = spkdata['spkindlst']
    spktlst          = spkdata['spktlst']
    num_cell_spks    = spkdata['num_cell_spks']
    pop_active_cells = spkdata['pop_active_cells']
    tmin             = spkdata['tmin']
    tmax             = spkdata['tmax']
    fraction_active  = { pop_name: float(len(pop_active_cells[pop_name])) / float(pop_num_cells[pop_name]) for pop_name in include }
    
    time_range = [tmin, tmax]

    # Calculate spike histogram if requested
    if spike_hist:
        all_spkts = np.concatenate([np.concatenate(lst, axis=0) for lst in spktlst])
        sphist_y, bin_edges = np.histogram(all_spkts, bins = np.arange(time_range[0], time_range[1], spike_hist_bin))
        sphist_x = bin_edges[:-1]+(spike_hist_bin / 2)

    maxN = 0
    minN = N

    avg_rates = {}
    tsecs = ((time_range[1]-time_range[0]) / 1e3)
    for i,pop_name in enumerate(spkpoplst):
        pop_num = len(pop_active_cells[pop_name])
        maxN = max(maxN, max(pop_active_cells[pop_name]))
        minN = min(minN, min(pop_active_cells[pop_name]))
        if pop_num > 0:
            if num_cell_spks[pop_name] == 0:
                avg_rates[pop_name] = 0
            else:
                avg_rates[pop_name] = ((num_cell_spks[pop_name] / pop_num) / tsecs)
        
    
    pop_colors = { pop_name: dflt_colors[ipop%len(raster_colors)] for ipop, pop_name in enumerate(spkpoplst) }

    pop_spk_dict = { pop_name: (pop_spkinds, pop_spkts) for (pop_name, pop_spkinds, pop_spkts) in zip(spkpoplst, spkindlst, spktlst) }

    if spike_hist is None:
        fig, axes = plt.subplots(nrows=len(spkpoplst), sharex=True, figsize=fig_options.figSize)
    elif spike_hist == 'subplot':
        fig, axes = plt.subplots(nrows=len(spkpoplst)+1, sharex=True, figsize=fig_options.figSize,
                                 gridspec_kw={'height_ratios': [1]*len(spkpoplst) + [2]})
    fig.suptitle ('CA1 Spike Raster', fontsize=fig_options.fontSize)

    sctplots = []
    
    for i, pop_name in enumerate(spkpoplst):

        if pop_name not in pop_spk_dict:
            continue
        
        pop_spkinds, pop_spkts = pop_spk_dict[pop_name]

        if max_spikes is not None:
            if int(max_spikes) < len(pop_spkinds):
               logger.info(('  Displaying only randomly sampled %i out of %i spikes for population %s' % (max_spikes, len(pop_spkts), pop_name)))
               sample_inds = np.random.randint(0, len(pop_spkinds)-1, size=int(max_spikes))
               pop_spkts   = pop_spkts[sample_inds]
               pop_spkinds = pop_spkinds[sample_inds]

        sct = None
        if len(pop_spkinds) > 0:
            sct = axes[i].scatter(pop_spkts, pop_spkinds, s=10, linewidths=fig_options.lw, marker=marker, c=pop_colors[pop_name], alpha=0.5, label=pop_name)
        axes[i].spines["top"].set_visible(False)
        axes[i].spines["bottom"].set_visible(False)
        axes[i].spines["left"].set_visible(False)
        axes[i].spines["right"].set_visible(False)
        sctplots.append(sct)

        N = pop_num_cells[pop_name]
        S = pop_start_inds[pop_name]
        axes[i].set_ylim(S, S+N-1)
        
    lgd_info = [(100. * fraction_active.get(pop_name, 0.), avg_rates.get(pop_name, 0.))
                for pop_name in spkpoplst ]
            
    # set raster plot y tick labels to the middle of the index range for each population
    for pop_name, a in zip_longest(spkpoplst, fig.axes[:-1]):
        if pop_name not in pop_active_cells:
            continue
        if len(pop_active_cells[pop_name]) > 0:
            maxN = max(pop_active_cells[pop_name])
            minN = min(pop_active_cells[pop_name])
            loc = pop_start_inds[pop_name] + 0.5 * (maxN - minN)
            yaxis = a.get_yaxis()
            yaxis.set_ticks([loc])
            yaxis.set_ticklabels([pop_name])
            yaxis.set_tick_params(length=0)
            a.get_xaxis().set_tick_params(length=0)
        
    # Plot spike histogram
    pch = interpolate.pchip(sphist_x, sphist_y)
    res_npts = int((sphist_x.max() - sphist_x.min()))
    sphist_x_res = np.linspace(sphist_x.min(), sphist_x.max(), res_npts, endpoint=True)
    sphist_y_res = pch(sphist_x_res)

    if spike_hist == 'overlay':
        ax2 = axes[-1].twinx()
        ax2.plot (sphist_x_res, sphist_y_res, linewidth=0.5)
        ax2.set_ylabel('Spike count', fontsize=fig_options.fontSize) # add yaxis label in opposite side
        ax2.set_xlim(time_range)
    elif spike_hist == 'subplot':
        ax2=axes[-1]
        ax2.plot (sphist_x_res, sphist_y_res, linewidth=1.0)
        ax2.set_xlabel('Time (ms)', fontsize=fig_options.fontSize)
        ax2.set_ylabel('Spikes', fontsize=fig_options.fontSize)
        ax2.set_xlim(time_range)
        
#    locator=MaxNLocator(prune='both', nbins=10)
#    ax2.xaxis.set_major_locator(locator)
    
    if labels == 'legend':
        # Shrink axes by 15%
        for ax in axes:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        if pop_rates:
            lgd_labels = [ '%s (%.02f%% active; %.3g Hz)' % (pop_name, info[0], info[1]) for pop_name, info in zip_longest(spkpoplst, lgd_info) ]
        else:
            lgd_labels = [ '%s (%.02f%% active)' % (pop_name, info[0]) for pop_name, info in zip_longest(spkpoplst, lgd_info) ]
        # Add legend
        lgd = fig.legend(sctplots, lgd_labels, loc = 'center right', 
                         fontsize='small', scatterpoints=1, markerscale=5.,
                         bbox_to_anchor=(1.002, 0.5), bbox_transform=plt.gcf().transFigure)
        fig.artists.append(lgd)
       
    elif labels == 'overlay':
        if pop_rates:
            lgd_labels = [ '%s (%.02f%% active; %.3g Hz)' % (pop_name, info[0], info[1]) for pop_name, info in zip_longest(spkpoplst, lgd_info) ]
        else:
            lgd_labels = [ '%s (%.02f%% active)' % (pop_name, info[0]) for pop_name, info in zip_longest(spkpoplst, lgd_info) ]
        for i, lgd_label in enumerate(lgd_labels):
            at = AnchoredText(pop_name + ' ' + lgd_label,
                              loc='upper right', borderpad=0.01, prop=dict(size=fig_options.fontSize))
            axes[i].add_artist(at)
        max_label_len = max([len(l) for l in lgd_labels])
        
    elif labels == 'yticks':
        for pop_name, info, a in zip_longest(spkpoplst, lgd_info, fig.axes[:-1]):
            if pop_rates:
                label = '%.02f%%\n%.2g Hz' % (info[0], info[1])
            else:
                label = '%.02f%%\n' % (info[0])

            maxN = max(pop_active_cells[pop_name])
            minN = min(pop_active_cells[pop_name])
            loc = pop_start_inds[pop_name] + 0.5 * (maxN - minN)
            a.set_yticks([loc, loc])
            a.set_yticklabels([pop_name, label])
            yticklabels = a.get_yticklabels()
            # Create offset transform in x direction
            dx = -66/72.; dy = 0/72. 
            offset = mpl.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
            # apply offset transform to labels.
            yticklabels[0].set_transform(yticklabels[0].get_transform() + offset)
            dx = -55/72.; dy = 0/72. 
            offset = mpl.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
            yticklabels[1].set_ha('left')    
            yticklabels[1].set_transform(yticklabels[1].get_transform() + offset)

            
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    
    # save figure
    if fig_options.saveFig:
       if isinstance(fig_options.saveFig, str):
           filename = fig_options.saveFig
       else:
           filename = f'{namespace_id} raster.{fig_options.figFormat}'
           plt.savefig(filename)
                
    # show fig 
    if fig_options.showFig:
        show_figure()
    
    return fig


def plot_lfp(input_path, config_path=None, time_range = None, compute_psd=False, window_size=4096, frequency_range=(0, 400.), overlap=0.9, bandpass_filter=False, dt=None, **kwargs):
    '''
    Line plot of LFP state variable (default: v). Returns figure handle.

    config: path to model configuration file
    input_path: file with LFP trace data
    time_range ([start:stop]): Time range of spikes shown; if None shows all (default: None)
    '''
    fig_options = copy.copy(default_fig_options)
    fig_options.update(kwargs)

    env = None
    if config_path is not None:
        env = Env(config_file=config_path)
    
    nrows = 1
    if env is not None:
        nrows = len(env.LFP_config)
    ncols = 1
    psd_col = 1
    if compute_psd:
        ncols += 1

    gs  = gridspec.GridSpec(nrows, ncols, width_ratios=[3,1] if ncols > 1 else [1])
    fig = plt.figure(figsize=fig_options.figSize)
    if env is None:

        lfp_array = np.loadtxt(input_path, dtype=np.dtype([("t", np.float32),
                                                           ("v", np.float32)]))

        if time_range is None:
            t = lfp_array['t']
            v = lfp_array['v']
        else:
            tlst = []
            vlst = []
            for (t,v) in zip(lfp_array['t'], lfp_array['v']):
                if time_range[0] <= t <= time_range[1]:
                    tlst.append(t)
                    vlst.append(v)
            t = np.asarray(tlst)
            v = np.asarray(vlst)

        if dt is None:
            raise RuntimeError("plot_lfp: dt must be provided when config_path is None")
        Fs = 1000. / dt

        if compute_psd:
            psd, freqs, peak_index = signal_psd(v, frequency_range=frequency_range, Fs=Fs, window_size=window_size, overlap=overlap)

        filtered_v = None
        if bandpass_filter:
            filtered_v = apply_filter(v, butter_bandpass(max(bandpass_filter[0], 1.0), bandpass_filter[1], Fs, order=2))

        iplot=0
        ax = plt.subplot(gs[iplot,0])
        ax.set_title('LFP', fontsize=fig_options.fontSize)
        ax.plot(t, v, linewidth=fig_options.lw)
        ax.set_xlabel('Time (ms)', fontsize=fig_options.fontSize)
        ax.set_ylabel('Field Potential (mV)', fontsize=fig_options.fontSize)
        
        if bandpass_filter:
            if filtered_v is not None:
                ax.plot(t, filtered_v, label='%s (filtered)' % lfp_label,
                        color='red', linewidth=fig_options.lw)
        if compute_psd:
            ax = plt.subplot(gs[iplot,psd_col])
            ax.plot(freqs, psd, linewidth=fig_options.lw)
            ax.set_xlabel('Frequency (Hz)', fontsize=fig_options.fontSize)
            ax.set_ylabel('Power Spectral Density (dB/Hz)', fontsize=fig_options.fontSize)
            ax.set_title('PSD (peak: %.3g Hz)' % (freqs[peak_index]), fontsize=fig_options.fontSize)
            
        if fig_options.saveFig:
            if isinstance(fig_options.saveFig, str):
                filename = fig_options.saveFig
            else:
                filename = f'CA1 LFP.{fig_options.figFormat}'
                plt.savefig(filename)
                
        # show fig
        if fig_options.showFig:
            show_figure()

    else:
        for iplot, (lfp_label, lfp_config_dict) in enumerate(viewitems(env.LFP_config)):
            namespace_id = "Local Field Potential %s" % str(lfp_label)
            import h5py
            infile = h5py.File(input_path)

            logger.info('plot_lfp: reading data for %s...' % namespace_id)
            if time_range is None:
                t = infile[namespace_id]['t']
                v = infile[namespace_id]['v']
            else:
                tlst = []
                vlst = []
                for (t,v) in zip(infile[namespace_id]['t'], infile[namespace_id]['v']):
                    if time_range[0] <= t <= time_range[1]:
                        tlst.append(t)
                        vlst.append(v)
                t = np.asarray(tlst)
                v = np.asarray(vlst)

            dt = lfp_config_dict['dt']
            Fs = 1000. / dt
        
            if compute_psd:
                freqs, psd, peak_index = signal_psd(v, Fs=Fs, frequency_range=frequency_range, window_size=window_size, overlap=overlap)

            filtered_v = None
            if bandpass_filter:
                filtered_v = apply_filter(v, butter_bandpass(max(bandpass_filter[0], 1.0), bandpass_filter[1], Fs, order=2))
                
            ax = plt.subplot(gs[iplot,0])
            ax.set_title('%s' % (namespace_id), fontsize=fig_options.fontSize)
            ax.plot(t, v, label=lfp_label, linewidth=fig_options.lw)
            ax.set_xlabel('Time (ms)', fontsize=fig_options.fontSize)
            ax.set_ylabel('Field Potential (mV)', fontsize=fig_options.fontSize)
            if bandpass_filter:
                if filtered_v is not None:
                    ax.plot(t, filtered_v, label='%s (filtered)' % lfp_label,
                            color='red', linewidth=fig_options.lw)
            if compute_psd:
                ax = plt.subplot(gs[iplot,psd_col])
                ax.plot(freqs, psd, linewidth=fig_options.lw)
                ax.set_xlabel('Frequency (Hz)', fontsize=fig_options.fontSize)
                ax.set_ylabel('Power Spectral Density (dB/Hz)', fontsize=fig_options.fontSize)
                ax.set_title('PSD (peak: %.3g Hz)' % (freqs[peak_index]), fontsize=fig_options.fontSize)

        # save figure
        if fig_options.saveFig:
            if isinstance(fig_options.saveFig, str):
                filename = fig_options.saveFig
            else:
                filename = namespace_id+'.%s' % fig_options.figFormat
                plt.savefig(filename)
                
        # show fig
        if fig_options.showFig:
            show_figure()

    return fig
