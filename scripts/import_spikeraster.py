import os, sys
import click
from ca1 import env, utils, io_utils

@click.command()
@click.option("--celltype-path", '-c', required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--input-path", '-i', required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--output-path", '-o',  type=click.Path(exists=False, file_okay=True, dir_okay=False))
@click.option("--output-namespace", '-n',  type=str, default="Spike Data")
@click.option("--output-npy", type=bool, default=False, is_flag=True)
@click.option("--verbose", "-v", type=bool, default=True, is_flag=True)
@click.option("--progress", "-p", type=bool, default=False, is_flag=True)
def main(celltype_path, input_path, output_path, output_namespace, output_npy, verbose, progress):

    utils.config_logging(verbose)
    io_utils.import_spikeraster(celltype_path, input_path, output_path, namespace=output_namespace, output_npy=output_npy, progress=progress)


if __name__ == '__main__':
    main(args=sys.argv[(utils.list_find(lambda x: os.path.basename(x) == os.path.basename(__file__), sys.argv)+1):])
