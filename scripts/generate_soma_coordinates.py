##
## Generate soma coordinates within layer-specific volume.
##

import os, sys, os.path, itertools, random, pickle, logging, click, gc
import math
from mpi4py import MPI
import h5py
import numpy as np
from neuroh5.io import append_cell_attributes, read_population_ranges
import rbf
from rbf.pde.geometry import contains
from rbf.pde.nodes import min_energy_nodes
from ca1.env import Env
from neural_geometry.alphavol import alpha_shape
from neural_geometry.geometry import make_uvl_distance, make_alpha_shape, load_alpha_shape, save_alpha_shape, get_total_extents, get_layer_extents, uvl_in_bounds
from ca1.CA1_volume import make_CA1_volume, CA1_volume
from ca1.utils import get_script_logger, config_logging, list_find, viewitems

script_name = os.path.basename(__file__)
logger = get_script_logger(script_name)

def mpi_excepthook(type, value, traceback):
    """

    :param type:
    :param value:
    :param traceback:
    :return:
    """
    sys_excepthook(type, value, traceback)
    if MPI.COMM_WORLD.size > 1:
        MPI.COMM_WORLD.Abort(1)


sys_excepthook = sys.excepthook
sys.excepthook = mpi_excepthook

def random_subset( iterator, K ):
    result = []
    N = 0

    for item in iterator:
        N += 1
        if len( result ) < K:
            result.append( item )
        else:
            s = int(random.random() * N)
            if s < K:
                result[ s ] = item

    return result


@click.command()
@click.option("--config", required=True, type=str)
@click.option("--config-prefix", required=False, type=click.Path(exists=True, file_okay=False, dir_okay=True), default="config")
@click.option("--types-path", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--geometry-path", required=False, type=click.Path(exists=False, file_okay=True, dir_okay=False))
@click.option("--output-path", required=True, type=click.Path(exists=False, file_okay=True, dir_okay=False))
@click.option("--output-namespace", type=str, default='Generated Coordinates')
@click.option("--populations", '-i', type=str, multiple=True)
@click.option("--resolution", type=(int,int,int), default=(3,3,3))
@click.option("--alpha-radius", type=float, default=1500.)
@click.option("--nodeiter", type=int, default=10)
@click.option("--optiter", type=int, default=200)
@click.option("--dispersion-delta", type=float, default=0.1)
@click.option("--snap-delta", type=float, default=0.01)
@click.option("--io-size", type=int, default=-1)
@click.option("--chunk-size", type=int, default=1000)
@click.option("--value-chunk-size", type=int, default=1000)
@click.option("--verbose", '-v', type=bool, default=False, is_flag=True)
def main(config, config_prefix, types_path, geometry_path, output_path, output_namespace, populations, resolution, alpha_radius, nodeiter, optiter, dispersion_delta, snap_delta, io_size, chunk_size, value_chunk_size, verbose):

    config_logging(verbose)
    logger = get_script_logger(script_name)

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    
    if io_size == -1:
        io_size = comm.size
    if rank == 0:
        logger.info('%i ranks have been allocated' % comm.size)


    if rank==0:
        if not os.path.isfile(output_path):
            input_file  = h5py.File(types_path,'r')
            output_file = h5py.File(output_path,'w')
            input_file.copy('/H5Types',output_file)
            input_file.close()
            output_file.close()
    comm.barrier()

    env = Env(comm=comm, config_file=config, config_prefix=config_prefix)

    random_seed = int(env.model_config['Random Seeds']['Soma Locations'])
    random.seed(random_seed)
    
    layer_extents = env.geometry['Parametric Surface']['Layer Extents']
    rotate = env.geometry['Parametric Surface']['Rotation']

    (extent_u, extent_v, extent_l) = get_total_extents(layer_extents)
    vol = make_CA1_volume(extent_u, extent_v, extent_l,
                          rotate=rotate, resolution=resolution)

    layer_alpha_shapes = {}
    layer_alpha_shape_path = 'Layer Alpha Shape/%d/%d/%d' % resolution
    if rank == 0:
        for layer, extents in viewitems(layer_extents):
            gc.collect()
            has_layer_alpha_shape = False
            if geometry_path:
                this_layer_alpha_shape_path = '%s/%s' % (layer_alpha_shape_path, layer)
                this_layer_alpha_shape = load_alpha_shape(geometry_path, this_layer_alpha_shape_path)
                layer_alpha_shapes[layer] = this_layer_alpha_shape
                if this_layer_alpha_shape is not None:
                    has_layer_alpha_shape = True
            if not has_layer_alpha_shape:
                logger.info("Constructing alpha shape for layers %s: extents: %s..." % (layer, str(extents)))
                #(extent_u, extent_v, extent_l) = get_total_extents(layer_extents)
                (extent_u, extent_v, extent_l) = get_layer_extents(layer_extents, layer)
                layer_vol = make_CA1_volume(extent_u, extent_v, extent_l,
                                      rotate=rotate, resolution=resolution)

                this_layer_alpha_shape = make_alpha_shape(layer_vol, alpha_radius=alpha_radius)
                layer_alpha_shapes[layer] = this_layer_alpha_shape
                if geometry_path:
                    save_alpha_shape(geometry_path, this_layer_alpha_shape_path, this_layer_alpha_shape)

    comm.barrier()
    population_ranges = read_population_ranges(output_path, comm)[0]
    if len(populations) == 0:
        populations = sorted(population_ranges.keys())

    if rank == 0:
        color = 1
    else:
        color = 0

    ## comm0 includes only rank 0
    comm0 = comm.Split(color, 0)

    for population in populations:

        (population_start, population_count) = population_ranges[population]

        pop_layers     = env.geometry['Cell Distribution'][population]
        pop_constraint = None
        if 'Cell Constraints' in env.geometry:
            if population in env.geometry['Cell Constraints']:
                pop_constraint = env.geometry['Cell Constraints'][population]
        if rank == 0:
            logger.info("Population %s: layer distribution is %s" % (population, str(pop_layers)))
            
        pop_layer_count = 0
        for layer, count in viewitems(pop_layers):
            pop_layer_count += count
        assert(population_count == pop_layer_count)


        xyz_coords = None
        xyz_coords_interp = None
        uvl_coords_interp = None
        if rank == 0:

            xyz_coords_lst = []
            xyz_coords_interp_lst = []
            uvl_coords_interp_lst = []
            for layer, count in viewitems(pop_layers):
                if count <= 0:
                    continue
                
                alpha = layer_alpha_shapes[layer]

                vert = alpha.points
                smp  = np.asarray(alpha.bounds, dtype=np.int64)
                for (vvi,vv) in enumerate(vert):
                    for (vi,v) in enumerate(vv):
                        if v <= 0.0: 
                            vert[vvi][vi] = 0.0

                N = int(count*2) # layer-specific number of nodes
                node_count = 0

                logger.info("Generating %i nodes..." % N)

                if verbose:
                    rbf_logger = logging.Logger.manager.loggerDict['rbf.pde.nodes']
                    rbf_logger.setLevel(logging.DEBUG)
                from rbf.pde.sampling import rejection_sampling
                while node_count < count:
                    # create N quasi-uniformly distributed nodes
                    def rho(x):
                        return np.ones(x.shape[0])
                    #nodes = rejection_sampling(N, rho, (vert, smp), start=0)
                    
                    out = min_energy_nodes(N,(vert,smp),iterations=nodeiter, 
                                           **{'dispersion_delta':dispersion_delta, 'snap_delta': snap_delta})
                    nodes = out[0]

                    # remove nodes with nan
                    nodes1 = nodes[~np.logical_or.reduce((np.isnan(nodes[:,0]), np.isnan(nodes[:,1]), np.isnan(nodes[:,2])))]
                    
                    # remove nodes outside of the domain
                    in_nodes = nodes[contains(nodes1, vert, smp)]
                    valid_idxs = None
                    if pop_constraint is not None:
                        if layer in pop_constraint:
                            valid_idxs = []
                            current_xyz = in_nodes.reshape(-1,3)
                            for i in range(len(current_xyz)):
                                if current_xyz[i][2] >= pop_constraint[layer][0] and current_xyz[i][2] <= pop_constraint[layer][1]:
                                    valid_idxs.append(i)
                            in_nodes = in_nodes[valid_idxs]
                    node_count = len(in_nodes)
                    N = int(1.5*N)
                    logger.info("%i interior nodes out of %i nodes generated" % (node_count, len(nodes)))

                xyz_coords_lst.append(in_nodes.reshape(-1,3))

            xyz_coords = np.concatenate(xyz_coords_lst)
            logger.info("Inverse interpolation of %i nodes..." % len(xyz_coords))
            uvl_coords_interp = vol.inverse(xyz_coords)
            xyz_coords_interp = vol(uvl_coords_interp[:,0],uvl_coords_interp[:,1],uvl_coords_interp[:,2],mesh=False).reshape(3,-1).T
            logger.info("Broadcasting generated nodes...")

            
        xyz_coords = comm.bcast(xyz_coords, root=0)
        xyz_coords_interp = comm.bcast(xyz_coords_interp, root=0)
        uvl_coords_interp = comm.bcast(uvl_coords_interp, root=0)

        coords = []
        coords_dict = {}
        xyz_error = np.asarray([0.0, 0.0, 0.0])

        if rank == 0:
            logger.info("Computing UVL coordinates...")

        for i in range(0,xyz_coords.shape[0]):

            coord_ind = i
            if i % size == rank:

                if uvl_in_bounds(uvl_coords_interp[coord_ind,:], layer_extents, pop_layers):
                    uvl_coords  = uvl_coords_interp[coord_ind,:].ravel()
                    xyz_coords1 = xyz_coords_interp[coord_ind,:].ravel()
                else:
                    uvl_coords = None
                    xyz_coords1 = None

                if uvl_coords is not None:

                    xyz_error   = np.add(xyz_error, np.abs(np.subtract(xyz_coords[coord_ind,:], xyz_coords1)))

                    logger.info('Rank %i: cell %i: %f %f %f' % (rank, i, uvl_coords[0], uvl_coords[1], uvl_coords[2]))

                    coords.append((xyz_coords1[0],xyz_coords1[1],xyz_coords1[2],
                                  uvl_coords[0],uvl_coords[1],uvl_coords[2]))
                                       
        
        total_xyz_error = np.zeros((3,))
        comm.Allreduce(xyz_error, total_xyz_error, op=MPI.SUM)

        coords_count = 0
        coords_count = np.sum(np.asarray(comm.allgather(len(coords))))

        if rank == 0:
            logger.info('Total %i coordinates generated' % coords_count)

        mean_xyz_error = np.asarray([(total_xyz_error[0] / coords_count), \
                                     (total_xyz_error[1] / coords_count), \
                                     (total_xyz_error[2] / coords_count)])

        
        if rank == 0:
            logger.info("mean XYZ error: %f %f %f " % (mean_xyz_error[0], mean_xyz_error[1], mean_xyz_error[2]))


        coords_lst = comm.gather(coords, root=0)
        if rank == 0:
            all_coords = []
            for sublist in coords_lst:
                for item in sublist:
                    all_coords.append(item)

            if coords_count < population_count:
                logger.warning("Generating additional %i coordinates " % (population_count - len(all_coords)))

                safety = 0.01
                delta = population_count - len(all_coords)
                for i in range(delta):
                    for layer, count in viewitems(pop_layers):
                        if count > 0:
                            min_extent = layer_extents[layer][0]
                            max_extent = layer_extents[layer][1]
                            coord_u = np.random.uniform(min_extent[0] + safety, max_extent[0] - safety)
                            coord_v = np.random.uniform(min_extent[1] + safety, max_extent[1] - safety)
                            coord_l = np.random.uniform(min_extent[2] + safety, max_extent[2] - safety)
                            xyz_coords = CA1_volume(coord_u, coord_v, coord_l, rotate=rotate).ravel()
                            all_coords.append((xyz_coords[0],xyz_coords[1],xyz_coords[2],\
                                              coord_u, coord_v, coord_l))

            sampled_coords = random_subset(all_coords, int(population_count))

            
            sampled_coords.sort(key=lambda coord: coord[3]) ## sort on U coordinate
            
            

            coords_dict = { population_start+i :  { 'X Coordinate': np.asarray([x_coord],dtype=np.float32),
                                    'Y Coordinate': np.asarray([y_coord],dtype=np.float32),
                                    'Z Coordinate': np.asarray([z_coord],dtype=np.float32),
                                    'U Coordinate': np.asarray([u_coord],dtype=np.float32),
                                    'V Coordinate': np.asarray([v_coord],dtype=np.float32),
                                    'L Coordinate': np.asarray([l_coord],dtype=np.float32) }
                            for (i,(x_coord,y_coord,z_coord,u_coord,v_coord,l_coord)) in enumerate(sampled_coords) }

            append_cell_attributes(output_path, population, coords_dict,
                                    namespace=output_namespace,
                                    io_size=io_size, chunk_size=chunk_size,
                                    value_chunk_size=value_chunk_size,comm=comm0)

        comm.barrier()

    comm0.Free()

if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda x: os.path.basename(x) == os.path.basename(__file__), sys.argv)+1):])
