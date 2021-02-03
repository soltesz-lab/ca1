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
from ca1.utils import get_module_logger, Struct, viewitems
from ca1.io_utils import get_h5py_attr, set_h5py_attr

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



def plot_coords_in_volume(populations, coords_path, coords_namespace, config, scale=25., subvol=False, verbose=False):
    
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
    for population in populations:
        coords = read_cell_attributes(coords_path, population, namespace=coords_namespace)

        count = 0
        for (k,v) in coords:
            count += 1
            xcoords.append(v['X Coordinate'][0])
            ycoords.append(v['Y Coordinate'][0])
            zcoords.append(v['Z Coordinate'][0])

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

    from mayavi import mlab
    
    logger.info('Plotting coordinates...')

    mlab.points3d(*pts.T, color=(1, 1, 0), scale_factor=scale)

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
        subvol.mplot_surface(color=(0, 0.4, 0), opacity=0.33)
    else:
        vol.mplot_surface(color=(0, 1, 0), opacity=0.33)
    
    mlab.show()
