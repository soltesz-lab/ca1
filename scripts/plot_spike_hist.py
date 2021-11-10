
import gc, os, sys
import click
import ca1
from ca1 import plot, utils

script_name = os.path.basename(__file__)


@click.command()
@click.option("--spike-events-path", '-p', required=True, type=click.Path())
@click.option("--spike-events-namespace", '-n', type=str, default='Spike Events')
@click.option("--populations", '-i', type=str, multiple=True)
@click.option("--include-artificial/--exclude-artificial", is_flag=True)
@click.option("--bin-size", type=float, default=5.0)
@click.option("--smooth", type=int)
@click.option("--t-variable", type=str, default='t')
@click.option("--t-max", type=float)
@click.option("--t-min", type=float)
@click.option("--quantity", type=str, default='rate')
@click.option("--font-size", type=float, default=14)
@click.option("--graph-type", type=str, default='bar')
@click.option("--overlay", is_flag=True)
@click.option("--save-format", type=str, default='png')
@click.option("--progress",  is_flag=True)
@click.option("--verbose", "-v", is_flag=True)
def main(spike_events_path, spike_events_namespace, populations, include_artificial, bin_size, smooth, t_variable, t_max, t_min, quantity, font_size, graph_type, overlay, save_format, progress, verbose):

    utils.config_logging(verbose)

    if t_max is None:
        time_range = None
    else:
        if t_min is None:
            time_range = [0.0, t_max]
        else:
            time_range = [t_min, t_max]

    if not populations:
        populations = ['eachPop']
        
    plot.plot_spike_histogram (spike_events_path, spike_events_namespace, include=populations, time_variable=t_variable,
                               time_range=time_range, pop_rates=True, bin_size=bin_size,
                               smooth=smooth, quantity=quantity, fontSize=font_size, overlay=overlay, graph_type=graph_type,
                               progress=progress, include_artificial=include_artificial, saveFig=True, figFormat=save_format)
    

if __name__ == '__main__':
    main(args=sys.argv[(utils.list_find(lambda x: os.path.basename(x) == script_name, sys.argv)+1):])
