
import os, os.path
try:
    from mpi4py import MPI  # Must come before importing NEURON
except Exception:
    pass
from ca1.utils import *
from neuron import h
from scipy import interpolate


def load_cell_template(env, pop_name, bcast_template=False):
    """
    :param pop_name: str
    """
    if pop_name in env.template_dict:
        return env.template_dict[pop_name]
    rank = env.comm.Get_rank()
    if not (pop_name in env.celltypes):
        raise KeyError('load_cell_templates: unrecognized cell population: %s' % pop_name)
    
    template_name = env.celltypes[pop_name]['template']
    if 'template file' in env.celltypes[pop_name]:
        template_file = env.celltypes[pop_name]['template file']
    else:
        template_file = None
    if not hasattr(h, template_name):
        find_template(env, template_name, template_file=template_file, path=env.template_paths, bcast_template=bcast_template)
    assert (hasattr(h, template_name))
    template_class = getattr(h, template_name)
    env.template_dict[pop_name] = template_class
    return template_class



def find_template(env, template_name, path=['templates'], template_file=None, bcast_template=False, root=0):
    """
    Finds and loads a template located in a directory within the given path list.
    :param env: :class:'Env'
    :param template_name: str; name of hoc template
    :param path: list of str; directories to look for hoc template
    :param template_file: str; file_name containing definition of hoc template
    :param root: int; MPI.COMM_WORLD.rank
    """
    if env.comm is None:
        bcast_template = False
    rank = env.comm.rank if env.comm is not None else 0
    found = False
    template_path = ''
    if template_file is None:
        template_file = '%s.hoc' % template_name
    if bcast_template:
        env.comm.barrier()
    if (env.comm is None) or (not bcast_template) or (bcast_template and (rank == root)):
        for template_dir in path:
            if template_file is None:
                template_path = '%s/%s.hoc' % (template_dir, template_name)
            else:
                template_path = '%s/%s' % (template_dir, template_file)
            found = os.path.isfile(template_path)
            if found and (rank == 0):
                logger.info('Loaded %s from %s' % (template_name, template_path))
                break
    if bcast_template:
        found = env.comm.bcast(found, root=root)
        env.comm.barrier()
    if found:
        if bcast_template:
            template_path = env.comm.bcast(template_path, root=root)
            env.comm.barrier()
        h.load_file(template_path)
    else:
        raise Exception('find_template: template %s not found: file %s; path is %s' %
                        (template_name, template_file, str(path)))


def configure_hoc_env(env, bcast_template=False):
    """
    :param env: :class:'Env'
    """
    h.load_file("stdrun.hoc")
    h.load_file("loadbal.hoc")
    for template_dir in env.template_paths:
        path = "%s/rn.hoc" % template_dir
        if os.path.exists(path):
            h.load_file(path)
    h.cvode.use_fast_imem(1)
    h.cvode.cache_efficient(1)
    h('objref pc, nc, nil')
    h('strdef dataset_path')
    if hasattr(env, 'dataset_path'):
        h.dataset_path = env.dataset_path if env.dataset_path is not None else ""
    if env.use_coreneuron:
        from neuron import coreneuron
        coreneuron.enable = True
        coreneuron.verbose = 1 if env.verbose else 0
    h.pc = h.ParallelContext()
    h.pc.gid_clear()
    env.pc = h.pc
    h.dt = env.dt
    h.tstop = env.tstop
    env.t_vec = h.Vector()  # Spike time of all cells on this host
    env.id_vec = h.Vector()  # Ids of spike times on this host
    env.t_rec = h.Vector() # Timestamps of intracellular traces on this host
    if 'celsius' in env.globals:
        h.celsius = env.globals['celsius']
    ## more accurate integration of synaptic discontinuities
    if hasattr(h, 'nrn_netrec_state_adjust'):
        h.nrn_netrec_state_adjust = 1
    ## sparse parallel transfer
    if hasattr(h, 'nrn_sparse_partrans'):
        h.nrn_sparse_partrans = 1