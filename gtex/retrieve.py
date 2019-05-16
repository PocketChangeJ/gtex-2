#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## file: retrieve.py
## desc: Retrieve GTEx eQTL and NCBI SNP datasets.
## auth: TR

from dask.distributed import Client
from dask.distributed import get_client
from dask.distributed import Future
from dask.distributed import LocalCluster
from typing import Dict
import gzip
import logging
import os
import requests as req
import shutil
import sys
import tarfile

from . import globe

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


def _download(url: str, output: str) -> None:
    """
    Download a file at the given URL and save it to disk.

    arguments
        url:    file URL
        output: output filepath
    """

    try:
        response = req.get(url, stream=True)

        ## Throw exception for any HTTP errors
        response.raise_for_status()

        with open(output, 'wb') as fl:
            for chunk in response.iter_content(chunk_size=1024):
                fl.write(chunk)

    except Exception as e:
        _logger.error('Request exception occurred: %s', e)
        raise


def _decompress(input: str, output: str) -> None:
    """
    Decompress a gzipped file.

    arguments
        input:  input filepath to a gzipped file
        output: output filepath
    """

    with gzip.open(input, 'rb') as gzfl, open(output, 'wb') as outfl:
        shutil.copyfileobj(gzfl, outfl)


def _extract_tar(input: str, output: str) -> None:
    """
    Extract all files from a TAR archive.

    arguments
        input:  input filepath to a TAR archive
        output: output filepath where all archives members are extracted
    """

    with tarfile.open(input) as tar:
        tar.extractall(output)


def download_gtex_eqtls(
    url: str = globe._url_eqtls,
    output: str = globe._fp_compressed_eqtls,
    force: bool = False
) -> None:
    """
    Retrieve GTEx eQTL datasets.

    arguments
        url:    optional GTEx eQTLs URL
        output: optional filepath where the compressed dataset would be stored
        force:  force data retrieval even if the dataset exists locally

    returns
        the output path string
    """

    if os.path.exists(output) and not force:
        _logger.warning('GTEx eQTL archive exists, skipping retrieval')
        return

    _logger.info('Retrieving GTEx eQTL archive')

    _download(url, output)


def download_gtex_lookup_table(
    url: str = globe._url_table,
    output: str = globe._fp_compressed_lookup_table,
    force: bool = False
) -> None:
    """
    Retrieve the GTEx eQTL lookup table.

    arguments
        url:    optional GTEx lookup table URL
        output: optional filepath where the compressed dataset would be stored
        force:  force data retrieval even if the dataset exists locally
    """

    if os.path.exists(output) and not force:
        _logger.warning('GTEx eQTL lookup table exists, skipping retrieval')
        return

    _logger.info('Retrieving GTEx lookup table')

    _download(url, output)


def download_gtex_annotations(
    url: str = globe._url_annotations,
    output: str = globe._fp_annotations,
    force: bool = False
) -> None:
    """
    Retrieve the GTEx annotations.

    arguments
        url:    optional GTEx annotations URL
        output: optional filepath where the dataset would be stored
        force:  force data retrieval even if the dataset exists locally
    """

    if os.path.exists(output) and not force:
        _logger.warning('GTEx eQTL annotations exist, skipping retrieval')
        return

    _logger.info('Retrieving GTEx annotations')

    _download(url, output)


def download_dbsnp_merge_table(
    url: str = globe._url_dbsnp150,
    output: str = globe._fp_compressed_dbsnp_table,
    force: bool = False
) -> None:
    """
    Retrieve the dbSNP merge table. Retrieves v150 by default.

    arguments
        url:    optional dbSNP merge table URL
        output: optional filepath where the compressed dataset would be stored
        force:  force data retrieval even if the dataset exists locally
    """

    if os.path.exists(output) and not force:
        _logger.warning('dbSNP merge table exists, skipping retrieval')
        return

    _logger.info('Retrieving NCBI dbSNP merge table')

    _download(url, output)


def decompress_gtex_lookup_table(
    input: str = globe._fp_compressed_lookup_table,
    output: str = globe._fp_lookup_table,
    **kwargs
) -> None:
    """
    Decompress the GTEx lookup table.

    arguments
        input:  optional filepath pointing to the compressed lookup table
        output: optional filepath where decompressed output is stored
        kwargs: the kwargs is just to trick dask for dependency tracking
    """

    _logger.info('Decompressing GTEx lookup table')

    _decompress(input, output)


def decompress_dbsnp_merge_table(
    input: str = globe._fp_compressed_dbsnp_table,
    output: str = globe._fp_dbsnp_table,
    **kwargs
) -> None:
    """
    Decompress the NCBI dbSNP merge table.

    arguments
        input:  optional filepath pointing to the compressed merge table
        output: optional filepath where decompressed output is stored
        kwargs: the kwargs is just to trick dask for dependency tracking
    """

    _logger.info('Decompressing NCBI dbSNP merge table')

    _decompress(input, output)


def extract_gtex_eqtls(
    input: str = globe._fp_compressed_eqtls,
    output: str = globe._dir_data_raw,
    **kwargs
) -> None:
    """
    Extract GTEx eQTLs to the given folder.

    arguments
        input:  optional filepath pointing to the eQTL archive
        output: optional filepath where extracted output is stored
        kwargs: the kwargs is just to trick dask for dependency tracking
    """

    _logger.info('Extracting GTEx eQTLs')

    _extract_tar(input, output)


def run_retrieve_step(client: Client = None) -> Dict[str, Future]:
    """
    Uses dask to run the data retrieval step of the ETL pipeline.

    arguments
        client: a dask Client object, if client is none then the function assumes
                it has been submitted to a worker node and will attempt to retrieve
                the client instance using get_client.

    returns
        a dict containing futures for each of the retrieved datasets
    """

    if not client:
        client = get_client()

    ## Download concurrently
    compressed_eqtls = client.submit(download_gtex_eqtls)
    compressed_annotations = client.submit(download_gtex_annotations)
    compressed_lookup_table = client.submit(download_gtex_lookup_table)
    compressed_merge_table = client.submit(download_dbsnp_merge_table)

    ## Decompress/extract concurrently
    eqtls = client.submit(extract_gtex_eqtls, depends=compressed_eqtls)
    lookup_table = client.submit(
        decompress_gtex_lookup_table, depends=compressed_lookup_table
    )
    merge_table = client.submit(
        decompress_dbsnp_merge_table, depends=compressed_merge_table
    )

    ## Wait for eQTL extraction to finish
    eqtls.result()

    tissue_eqtls = []

    ## Iterate through the directory and decompress/move eQTLs for each tissue
    for fl in os.listdir(globe._dir_eqtls):

        ## The directory contains two types of datasets: eQTLs and eGenes. We're only
        ## interested in the former.
        if 'signif_variant_gene_pairs.txt.gz' not in fl:
            continue

        input = os.path.join(globe._dir_eqtls, fl)
        output = os.path.join(globe._dir_data_raw, os.path.splitext(fl)[0])

        tissue_eqtls.append(client.submit(_decompress, input, output))

    ## Return the futures
    return {
        'eqtls': tissue_eqtls,
        'annotations': compressed_annotations,
        'lookup_table': lookup_table,
        'merge_table': merge_table
    }


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

    _initialize_logging(args.verbose)

    _logger.info('Starting cluster for retrieval step')

    ## Create a local cluster
    client = Client(LocalCluster(
        n_workers = 12,
        processes=True
    ))

    client.run(_initialize_logging, args.verbose)

    ## Run the data retrieval step
    futures = run_retrieve_step(client)

    ## Wait on the results
    client.gather(futures)

    _logger.info('Finished data retrieval')

    client.close()

