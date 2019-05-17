#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## file: process.py
## desc: Process GTEx eQTL SNP datasets.
## auth: TR

from dask.distributed import Client
from dask.distributed import get_client
from dask.distributed import Future
from dask.distributed import LocalCluster
from dask.distributed import wait
from pathlib import Path
from typing import Dict, List
import dask.dataframe as dd
import logging
import os
import numpy as np
import pandas as pd
import re
import sys

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


def parse_annotations(fp: str = globe._fp_annotations) -> pd.DataFrame:
    """
    Parse the GTEx sample annotations file. We're mainly interesting in getting tissue
    names/groups.

    arguments
        fp: optional filepath to the annotations

    returns
        a dataframe of tissue names and groups
    """

    df = pd.read_csv(fp, sep='\t')

    ## These are the tissue groups and names respectively
    df = df[['SMTS', 'SMTSD']]

    ## Rename columns
    df = df.rename(columns={'SMTS': 'tissue_group', 'SMTSD': 'tissue_name'})

    ## Drop duplicate tissues
    df = df.drop_duplicates(subset='tissue_name')

    ## Add columns that can be used for matching data filenames -> tissue
    df['tissue_str'] = df.tissue_name.str.replace('\(|\)| - ', ' ')
    df['tissue_str'] = df.tissue_str.str.replace('\s+', '.*')

    ## Drop tissue group from the tissue name
    df['tissue_name'] = df.tissue_name.str.replace('\w+ - ', '')

    return df


#def parse_lookup_table(fp: str) -> Dict[str, str]:
def parse_lookup_table(fp: str = globe._fp_lookup_table) -> pd.DataFrame:
    """
    The GTEx lookup table is tab delimited and we're interested in the following
    columns (indexes):

    (0) (2)        (6)
    chr variant_id rsid_id_dbSNP147_GRCh37p13

    The variant_id is the unique GTEx identifier which is also found in the eQTL data
    files. We use it to associate each eQTL with a canonical reference SNP identifier
    (rsID).

    Since this is a larger file (~1.8GB) it's faster to parse it as a dask dataframe then
    join it into a single pandas DF.

    arguments
        fp: filepath to the GTEx lookup table

    returns
        a mapping of GTEx variant IDs to dbSNP identifiers (rsID)
    """

    df = dd.read_csv(fp, blocksize='400MB', sep='\t', dtype={'chr': str})

    ## Rename ugly ass column name
    df = df.rename(columns={'rs_id_dbSNP147_GRCh37p13': 'rsid'})

    ## Remove anything that doesn't have an rsID
    df = df[df.rsid != '.']

    ## Only keep the GTEx variant ID and the rsID
    df = df[['variant_id', 'rsid']]

    ## Compute and return the pandas DF to the client node
    df = df.compute()

    return dict(zip(df.variant_id, df.rsid))
    ## Convert to a dict
    #return dict(zip(df.variant_id, df.rsid))
    """
    """

    """
    table = []

    for df in pd.read_csv(fp, sep='\t', chunksize=4096):
        #df = dd.read_csv(fp, blocksize='250MB', sep='\t', dtype={'chr': str})

        ## Rename ugly ass column name
        df = df.rename(columns={'rs_id_dbSNP147_GRCh37p13': 'rsid'})

        ## Remove anything that doesn't have an rsID
        df = df[df.rsid != '.']

        ## Only keep the GTEx variant ID and the rsID
        df = df[['variant_id', 'rsid']]

        table.extend(zip(df.variant_id, df.rsid))
        ## Compute and return the pandas DF to the client node
        #df = df.compute()

    return dict([(v, r) for v, r in table])
    """

    #return df
    #table.extend(zip(df.variant_id, df.rsid))

    #return dict([(v, r) for v, r in table])


def parse_eqtls(fp: str) -> pd.DataFrame:
    """
    Each tissue-specific eQTL file is tab delimited and we're interested in the
    following columns:

    (0)        (1)     (6)
    variant_id gene_id p-value

    arguments
        fp: filepath to the GTEx eQTL dataset

    returns
        a data frame containing the parsed file. The frame only contains the
        variant_id and gene_id columns.
    """

    df = pd.read_csv(fp, sep='\t')

    ## Remove the Ensembl gene version
    df['gene_id'] = df.gene_id.map(lambda g: g.split('.')[0])

    ## Remove anything w/ p < 0.05 (I don't think there are any but just in case)
    df = df[df.pval_nominal < 0.05]

    ## Rename p-value column
    df = df.rename(columns={'pval_nominal': 'p'})

    ## Add the filename to the frame which is used for tissue annotation later on
    df['filename'] = Path(fp).name

    return df[['variant_id', 'gene_id', 'p', 'filename']]


def parse_merge_table(fp: str = globe._fp_dbsnp_table) -> pd.DataFrame:
    """
    Parse the NCBI dbSNP merge table. A description of the table can be found here:
    https://www.ncbi.nlm.nih.gov/projects/SNP/snp_db_table_description.cgi?t=RsMergeArch

    arguments
        fp: filepath to the merge table

    returns
        a dataframe of the merge table
    """

    df = dd.read_csv(
        fp,
        blocksize='250MB',
        sep='\t',
        assume_missing=True,
        header=None,
        names=[
            'high',
            'low',
            'build',
            'orien',
            'created',
            'updated',
            'current',
            'o2c',
            'comment'
        ]
    )

    # We really only need high and current fields for the mapping procedure
    df = df[['high', 'current']]

    ## Return a pandas DF to the client node
    return df.compute()


def annotate_tissue(annotations: pd.DataFrame, eqtls: pd.DataFrame) -> pd.DataFrame:
    """
    Attempt to associate a tissue and tissue group with the given eQTL dataset.

    """

    ## Get the filename from the given set of eQTLs
    filename = eqtls.filename.loc[0]
    _logger.info('filename: %s', filename)

    for annotation in annotations.itertuples():
        #_logger.info('tissue_str: %s', annotation.tissue_str)
        ## Try to match the the filename and annotation
        if re.search(annotation.tissue_str, filename, re.IGNORECASE):
            _logger.info('found tissue group')
            eqtls['tissue'] = annotation.tissue_name
            eqtls['tissue_group'] = annotation.tissue_group
            break

    else:
        _logger.warning('Could not discern tissue type for %s', filename)

        eqtls['tissue'] = 'unknown'
        eqtls['tissue_group'] = 'unknown'

    _logger.info('dropping filename column')
    ## The filename column is no longer needed
    eqtls = eqtls.drop(columns='filename')

    _logger.info('returning')
    return eqtls


def map_eqtl_rsids(lookup: Dict[str, str], eqtls: pd.DataFrame) -> pd.DataFrame:
    """
    Map GTEx variant IDs to canonical reference SNP identifiers (rsID).

    arguments
        lookup: the GTEx lookup table mapping
        eqtls:  an eQTL dataset

    returns
    """

    _logger.info('mapping eqtl rsIDs')
    ## Map to dbSNP references
    eqtls['variant_id'] = eqtls.variant_id.map(lambda v: lookup.get(v))
    _logger.info('done mapping eqtl rsIDs')

    print(eqtls.head())
    return eqtls.iloc[:100]
    ## Keep tab on simple stats associated with each eQTL dataset (e.g. # of total SNPs,
    ## of SNPs that didn't map to rsIDs, etc.)
    #statistics = {
    #    'tissue': eqtls.tissue.iloc[0],
    #    'tissue_group': eqtls.tissue_group.iloc[0],
    #    'no_rsid': len(eqtls[eqtls.rsid.isnull()].index),
    #    'have_rsid': len(eqtls[eqtls.rsid.notnull()].index),
    #}

    #return (eqtls, statistics)


def merge_snps(merge: pd.DataFrame, eqtls: pd.DataFrame) -> pd.DataFrame:
    """
    Update SNP identifiers, where possible, for the given eQTL dataset using
    the dbSNP merge table.

    arguments
        merge: dataframe containing the dbSNP merge table
        eqtls: dataframe containing partially processed eQTLs

    returns
        a dataframe with updated SNP IDs
    """

    ## Strip out the rs prefix from the SNP identifier and convert the ID to an integer
    eqtls['rsid'] = eqtls['rsid'].str.strip(to_strip='rs')
    eqtls['rsid'] = eqtls['rsid'].astype(np.int64)

    ## Join the dataframes on the old SNP identifiers (i.e. 'high' in the merge table)
    eqtls = eqtls.join(merge.set_index('high'), on='rsid', how='left')

    return eqtls
    #criterion = eqtls.current.notna()

    ### Update IDs if necessary
    #eqtls.loc[criterion, args.column] = eqtls.loc[criterion, 'current']

    #eqtls = eqtls.drop(columns='current')

    ### Add the rs prefix back in
    #eqtls[args.column] = 'rs' + eqtls[args.column].astype(np.int64).astype('str')


def finalize_eqtls(eqtls: pd.DataFrame) -> pd.DataFrame:
    """

    :param eqtls:
    :return:
    """

    ## Replace rsIDs with the updated ones from the merge table
    criterion = eqtls.current.notnull()

    eqtls.loc[criterion, 'rsid'] = eqtls.loc[criterion, 'current']

    eqtls = eqtls.drop(columns='current')
    eqtls = eqtls.dropna(axis=0, subset=['rsid'])

    ## Add the rs prefix back
    eqtls['rsid'] = 'rs' + eqtls['rsid'].astype(np.int64).astype(str)

    return eqtls


def save_eqtls(eqtls: List[pd.DataFrame], outdir: str = globe._dir_data_processed) -> None:
    """
    :param eqtls:
    :param outdir:
    :return:
    """

    ## Merge eQTLs into a single dataframe
    eqtls = pd.concat(eqtls)

    ## Group by tissue group and save
    for (group, df) in eqtls.groupby(['tissue_group']):

        ## Format group name for output
        group = f"{group.replace(' ', '-')}.tsv"

        ## Output filepath
        outpath = Path(outdir).joinpath(group)

        ## Save the output
        df.to_csv(
            outpath,
            sep='\t',
            index=False,
            columns=['rsid', 'gene_id', 'p', 'tissue', 'tissue_group']
        )


def calculate_eqtl_stats(eqtls: pd.DataFrame) -> Dict[str, str]:
    """
    Calculate simple summary stats for a single eQTL dataset.

    arguments
        eqtls: a dataframe contaning partially processed eQTLs

    returns
        a dict with some summary stats
    """

    return {
        'tissue': eqtls.tissue.iloc[0],
        'tissue_group': eqtls.tissue_group.iloc[0],
        'no_rsid': str(len(eqtls[eqtls.rsid.isnull()].index)),
        'have_rsid': str(len(eqtls[eqtls.rsid.notnull()].index)),
        'new_rsid': str(len(eqtls[eqtls.current.notnull()].index)),
    }


def run_processing_step(client: Client = None) -> Dict[str, Future]:
    """
    Uses dask to run the data processing step of the GTEx ETL pipeline.

    arguments
        client: a dask Client object, if client is none then the function assumes
                it has been submitted to a worker node and will attempt to retrieve
                the client instance using get_client.

    returns
        a dict containing futures for each of the retrieved datasets
    """

    if not client:
        client = get_client()

    ## Lookup and merge table parsing isn't submitted since we use a dask dataframe
    ## to parse them
    annotations = client.submit(parse_annotations)
    #_logger.info(annotations.result().head())
    lookup = parse_lookup_table()
    lookup = client.scatter(lookup)
    #_logger.info(list(lookup.items())[:5])
    #return
    merge = parse_merge_table()

    future_eqtls = []
    future_stats = []

    ## For each tissue-specific eQTL dataset...
    for fp in Path(globe._dir_data_raw).iterdir():

        ## Make sure it's actually an eQTL file
        if 'signif_variant_gene_pairs.txt' not in fp.name:
            continue

        ## Parse the eQTLs
        eqtls = client.submit(parse_eqtls, fp)

        ## Annotate the tissue
        annotated_eqtls = client.submit(annotate_tissue, annotations, eqtls)
        #_logger.info(annotated_eqtls.result().head())
        #future_eqtls.append(annotated_eqtls)

        #annotated_eqtls = client.scatter(annotated_eqtls)

        ## Map variant IDs -> rsIDs
        mapped_eqtls = client.submit(map_eqtl_rsids, lookup, annotated_eqtls)
        wait(mapped_eqtls)
        #_logger.info(mapped_eqtls.result().head())

        ### Update SNP identifiers to their latest versions
        #merged_snps = client.submit(merge_snps, merge, mapped_eqtls)

        ### Calculate summary stats for output
        #stats = client.submit(calculate_eqtl_stats, merged_snps)

        #final_eqtls = client.submit(finalize_eqtls, merged_snps)

        #future_eqtls.append(final_eqtls)
        break


    client.gather(future_eqtls)
    ## Now gather all the results, concat them into a single DF, group by tissue and save
    #save_eqtls(client.gather(future_eqtls))
    #saved_eqtls = client.submit(save)


    ## Download concurrently
    #compressed_eqtls = client.submit(download_gtex_eqtls)
    #compressed_annotations = client.submit(download_gtex_annotations)
    #compressed_lookup_table = client.submit(download_gtex_lookup_table)
    #compressed_merge_table = client.submit(download_dbsnp_merge_table)

    ### Decompress/extract concurrently
    #eqtls = client.submit(extract_gtex_eqtls, depends=compressed_eqtls)
    #lookup_table = client.submit(
    #    decompress_gtex_lookup_table, depends=compressed_lookup_table
    #)
    #merge_table = client.submit(
    #    decompress_dbsnp_merge_table, depends=compressed_merge_table
    #)

    ### Wait for eQTL extraction to finish
    #eqtls.result()

    #tissue_eqtls = []

    ### Iterate through the directory and decompress/move eQTLs for each tissue
    #for fl in os.listdir(globe._dir_eqtls):

    #    ## The directory contains two types of datasets: eQTLs and eGenes. We're only
    #    ## interested in the former.
    #    if 'signif_variant_gene_pairs.txt.gz' not in fl:
    #        continue

    #    input = os.path.join(globe._dir_eqtls, fl)
    #    output = os.path.join(globe._dir_data_raw, os.path.splitext(fl)[0])

    #    tissue_eqtls.append(client.submit(_decompress, input, output))

    ### Return the futures
    #return {
    #    'eqtls': tissue_eqtls,
    #    'annotations': compressed_annotations,
    #    'lookup_table': lookup_table,
    #    'merge_table': merge_table
    #}


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

    _logger.info('Starting cluster for processing step')

    ## Create a local cluster
    client = Client(LocalCluster(
        n_workers=10,
        processes=True
    ))

    client.run(_initialize_logging, args.verbose)

    run_processing_step(client)
    ## Run the data retrieval step
    #futures = run_retrieve_step(client)

    ## Wait on the results
    #client.gather(futures)

    #_logger.info('Finished data retrieval')

    client.close()

