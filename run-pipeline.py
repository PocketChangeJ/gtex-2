#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## file: run-pipeline.py
## desc: Run the GTEx eQTL processing pipeline.
## auth: TR

from dask.distributed import Client
from dask.distributed import LocalCluster
import logging
import sys

from gtex import process
from gtex import retrieve

_logger = logging.getLogger(__name__)


def _initialize_logging(verbose: bool) -> None:
    """
    Initialize logging across workers.

    arguments
        verbose: indicates if verbose logging is on/off
    """

    ## Add a console logger based on verbosity settings
    conlog = logging.StreamHandler()
    conlog.setLevel(logging.INFO if verbose else logging.ERROR)
    conlog.setFormatter(logging.Formatter('[%(levelname)-7s] %(message)s'))

    _logger.setLevel(logging.INFO if verbose else logging.ERROR)
    _logger.addHandler(conlog)


if __name__ == '__main__':
    from argparse import ArgumentParser

    usage = f'{sys.argv[0]} [options]'
    parser = ArgumentParser(usage=usage)

    parser.add_argument(
        '--verbose',
        action='store_true',
        dest='verbose',
        help='clutter your screen with output'
    )

    args = parser.parse_args()

    _initialize_logging(args.verbose)

    ## Create a local cluster
    client = Client(LocalCluster(
        n_workers=8,
        processes=True
    ))

    client.run(_initialize_logging, args.verbose)

    _logger.info('Starting data retrieval')

    ## Run the data retrieval step
    futures = retrieve.run_retrieve_step(client)

    ## Wait for retrieval and extraction to finish
    client.gather(futures)
    exit()

    _logger.info('Starting eQTL processing')

    ## Process and save eQTL datasets
    process.run_processing_step()

    client.close()

