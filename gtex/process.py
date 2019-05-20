#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## file: process.py
## desc: Process GTEx eQTL SNP datasets.
## auth: TR

from pathlib import Path
from typing import Dict, List
import logging
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
    df.loc[:, 'tissue_str'] = df.tissue_name.str.replace('\(|\)| - ', ' ')
    df.loc[:, 'tissue_str'] = df.tissue_str.str.replace('\s+', '.*')

    ## Drop tissue group from the tissue name
    df.loc[:, 'tissue_name'] = df.tissue_name.str.replace('\w+ - ', '')

    return df


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

    df = pd.read_csv(fp, sep='\t', dtype={'chr': str})

    ## Rename ugly ass column name
    df = df.rename(columns={'rs_id_dbSNP147_GRCh37p13': 'rsid'})

    ## Remove anything that doesn't have an rsID
    df = df[df.rsid != '.']

    return dict(zip(df.variant_id, df.rsid))


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
    df.loc[:, 'gene_id'] = df.gene_id.map(lambda g: g.split('.')[0])

    ## Remove anything w/ p < 0.05 (I don't think there are any but just in case)
    df = df[df.pval_nominal < 0.05]

    ## Rename p-value column
    df = df.rename(columns={'pval_nominal': 'p'})

    ## Add the filename to the frame which is used for tissue annotation later on
    df.loc[:, 'filename'] = Path(fp).name

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

    df = pd.read_csv(
        fp,
        sep='\t',
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

    return dict(zip(df.high, df.current))


def annotate_tissue(
        annotations: pd.DataFrame,
        eqtls: pd.DataFrame,
        filename: str
) -> pd.DataFrame:
    """
    Attempt to associate a tissue and tissue group with the given eQTL dataset.

    arguments
        annotations: GTEx tissue annotations
        eqtls:       eQTL dataset
        filename:    filename of the eQTL dataset which should contain tissue info

    returns
        annotated eQTLs
    """

    for annotation in annotations.itertuples():
        ## Try to match the the filename and annotation
        if re.search(annotation.tissue_str, filename, re.IGNORECASE):
            eqtls.loc[:, 'tissue'] = annotation.tissue_name
            eqtls.loc[:, 'tissue_group'] = annotation.tissue_group
            break

    else:
        _logger.warning('Could not discern tissue type for %s', filename)

        eqtls.loc[:, 'tissue'] = 'unknown'
        eqtls.loc[:, 'tissue_group'] = 'unknown'

    eqtls = eqtls.drop('filename', axis=1)

    return eqtls


def map_eqtl_rsids(lookup: Dict[str, str], eqtls: pd.DataFrame) -> pd.DataFrame:
    """
    Map GTEx variant IDs to canonical reference SNP identifiers (rsID).

    arguments
        lookup: the GTEx lookup table mapping
        eqtls:  an eQTL dataset

    returns
        eQTLs with their rsID mappings
    """

    ## Map to dbSNP references by joining on the variant_id
    eqtls.loc[:, 'rsid'] = eqtls.variant_id.map(lambda v: lookup.get(v))

    ## Set missing rsIDs (should be NaN values) to zero
    eqtls = eqtls.fillna(value={'rsid': 'rs0'})

    return eqtls


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
    eqtls.loc[:, 'rsid'] = eqtls['rsid'].str.strip(to_strip='rs')
    eqtls.loc[:, 'rsid'] = eqtls['rsid'].astype(np.int64)

    ## Join the dataframes on the old SNP identifiers (i.e. 'high' in the merge table)
    eqtls.loc[:, 'merged'] = eqtls.rsid.map(lambda v: v in merge)
    eqtls.loc[:, 'rsid'] = eqtls.rsid.map(lambda v: merge.get(v, v))

    return eqtls


def finalize_eqtls(eqtls: pd.DataFrame) -> pd.DataFrame:
    """
    Finalizes an eQTL dataset by 1) removing SNPs with missing rsIDs, 2) adding the rs
    prefix to each ID, and 3) removing columns that are no longer needed.

    arguments
        eqtls: eQTL dataset

    returns
        eQTLs
    """

    eqtls = eqtls[eqtls.rsid != 0]

    ## Add the rs prefix back
    eqtls.loc[:, 'rsid'] = 'rs' + eqtls['rsid'].astype(np.int64).astype(str)

    return eqtls[['rsid', 'gene_id', 'p', 'tissue', 'tissue_group']]


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
        'no_rsid': str(len(eqtls[eqtls.rsid == 0].index)),
        'have_rsid': str(len(eqtls[eqtls.rsid != 0].index)),
        'new_rsid': str(len(eqtls[eqtls.merged].index)),
    }


def save_eqtls(
        eqtls: List[pd.DataFrame],
        outdir: str = globe._dir_data_processed
) -> None:
    """
    Group eQTLs by tissue and save the grouped eQTLs to a file.

    arguments
        eqtls: a list of eQTL datasets
        outdir: an output directory path
    """

    ## Merge eQTLs into a single dataframe
    eqtls = pd.concat(eqtls)

    ## Group by tissue group and save
    for (group, df) in eqtls.groupby(['tissue_group']):

        ## Format group name for output
        group = f"{group.replace(' ', '-')}.tsv".lower()

        ## Output filepath
        outpath = Path(outdir).joinpath(group)

        ## Save the output
        df.to_csv(
            outpath,
            sep='\t',
            index=False,
            columns=['rsid', 'gene_id', 'p', 'tissue', 'tissue_group']
        )


def save_eqtl_stats(stats: Dict[str, str], output: str) -> None:
    """
    Save eQTL mapping summary stats to a file.

    arguments
        stats:  eQTL stats
        output: output filepath
    """

    pd.DataFrame(stats).to_csv(
        output,
        sep='\t',
        index=False,
    )


def run_processing_step() -> None:
    """
    Run the data processing step of the GTEx ETL pipeline.
    Processes each file GTEx file individually and saves the results per tissue group.
    """

    _logger.info('Parsing annotations, GTEx lookup table, and dbSNP merge table')

    annotations = parse_annotations()
    lookup = parse_lookup_table()
    merge = parse_merge_table()

    processed_eqtls = []
    eqtl_stats = []

    _logger.info('Parsing eQTL datasets')

    ## For each tissue-specific eQTL dataset...
    for fp in Path(globe._dir_data_raw).iterdir():

        ## Make sure it's actually an eQTL file
        if 'signif_variant_gene_pairs.txt' not in fp.name:
            continue

        ## Parse the eQTLs
        eqtls = parse_eqtls(fp)

        ## Annotate the tissue
        annotated_eqtls = annotate_tissue(annotations, eqtls, fp.name)

        ## Map variant IDs -> rsIDs
        mapped_eqtls = map_eqtl_rsids(lookup, annotated_eqtls)

        ## Update SNP identifiers to their latest versions
        merged_snps = merge_snps(merge, mapped_eqtls)

        eqtl_stats.append(calculate_eqtl_stats(merged_snps))

        ## Calculate summary stats for output
        final_eqtls = finalize_eqtls(merged_snps)

        processed_eqtls.append(final_eqtls)

    _logger.info('Saving processed datasets')

    save_eqtl_stats(eqtl_stats, globe._fp_eqtl_stats)
    save_eqtls(processed_eqtls)


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

    run_processing_step()

