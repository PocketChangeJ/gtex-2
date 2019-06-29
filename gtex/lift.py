#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## file: bedops.py
## desc: Python wrappers for common bedops tools and operations.

from typing import List
from dask.distributed import Client
from dask.distributed import Future
from dask.distributed import LocalCluster
from dask.distributed import get_client
from functools import partial
from glob import glob
from pathlib import Path
import logging
import numpy as np
import pandas as pd
import subprocess as sub
import sys
import tempfile as tf

import time

from . import globe
from . import log

logging.getLogger(__name__).addHandler(logging.NullHandler())


def _update_genome_build_coordinates(
    df: pd.DataFrame,
    unmapped: str,
    chain: str = globe._fp_hg38_chain,
    liftover: str = globe._exe_liftover
) -> pd.DataFrame:
    """
    Update SNP coordinates to the latest genome build (hg38).

    arguments
        df:       eQTL dataframe
        unmapped: output filepath used to store unmapped coords
        chain:    chain file needed by the liftOver exe
        liftover: optional path to the liftOver exe
    """

    tmp_input = tf.NamedTemporaryFile(delete=True)
    tmp_output = tf.NamedTemporaryFile(delete=True)

    ## Generate a UID for each row so we can rejoin metadata later on
    df = df.reset_index(drop=False)
    df = df.rename(columns={'index': 'uid'})

    ## Prior to overlapping we have to add an end column to the eQTL DF otherwise
    ## liftOver will shit the bed
    df['end'] = df.start + 1

    ## Since liftOver is really fucking picky about input, we have to remove all extra
    ## columns and get rid of the header when saving to the temp input file
    df[['chromosome', 'start', 'end', 'uid']].to_csv(
        tmp_input.name, sep='\t', index=False, header=False
    )

    ## Run liftover, usage: liftOver oldFile map.chain newFile unMapped
    try:
        sub.run(
            [liftover, tmp_input.name, chain, tmp_output.name, unmapped],
            stderr=sub.DEVNULL
        )

    except Exception as e:
        log._logger.error('There was a problem running liftOver as a subprocess: %s', e)
        raise

    ## Now we need to read in the lifted coordinates
    lifted_df = pd.read_csv(
        tmp_output.name,
        sep='\t',
        header=None,
        names=['chromosome', 'start', 'end', 'uid']
    )

    ## Rejoin lifted coordinates back to the original metadata, ignoring anything that
    ## wasn't lifted to the latest build
    joined = lifted_df.join(
        df.drop(columns=['chromosome', 'start', 'end']).set_index('uid'),
        on='uid',
        how='left'
    )

    return joined.drop(columns=['uid'])


def _lift_genome_build(input: str, output: str, unmapped: str) -> str:
    """

    :param input:
    :param output:
    :param unmapped:
    :return:
    """

    df = pd.read_csv(input, sep='\t')
    lifted = _update_genome_build_coordinates(df, unmapped)

    lifted.to_csv(output, sep='\t', index=False)

    return output


def run_lift_eqtls(
    indir: str = globe._dir_data_processed,
    outdir: str = globe._dir_data_lifted
) -> List[Future]:
    """

    :param indir:
    :return:
    """

    client = get_client()

    futures = []

    for fp in glob(Path(indir, '*.tsv').as_posix()):

        output = Path(outdir, Path(fp).name)
        unmapped = Path(outdir, Path(fp).with_suffix('.unmapped').name)

        future = client.submit(
            _lift_genome_build, fp, output.as_posix(), unmapped.as_posix()
        )

        futures.append(future)

    return futures


## This step in the GTEx pipeline can be run individually
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

    ## Create a local cluster
    client = Client(LocalCluster(
        n_workers=10,
        processes=True,
        local_dir='/var/tmp'
    ))

    ## Run the logging init function on each worker and register the callback so
    ## future workers also run the function
    init_logging_partial = partial(log.initialize_logging, verbose=args.verbose)
    client.register_worker_callbacks(setup=init_logging_partial)
    log.initialize_logging(verbose=args.verbose)

    futures = run_lift_eqtls()

    client.gather(futures)

