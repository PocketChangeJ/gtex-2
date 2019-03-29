#!/usr/bin/env python

## file: merge_snp_ids.py
## desc: Use the NCBI RsMergeArch table to merge old SNP identifiers into their 
##       current versions.
## auth: TR

from __future__ import print_function
from sys import argv
import logging
import numpy as np
import os
import pandas as pd
import re

LOG = logging.getLogger(__name__)

def find_eqtl_files(d):
    """
    Attempts to locate eQTL variant-gene pair files in the given directory.

    arguments
        d: some directory

    returns
        a list of eQTL filepaths
    """

    return [
        os.path.join(d, f) for f in os.listdir(d) if 'signif_variant_gene_pairs' in f
    ]

def parse_annotations_file(fp):
    """
    Parse the GTEx sample annotations file. We're mainly interesting in getting tissue
    names/groups.

    arguments
        fp: filepath to the annotations

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

def parse_lookup_file(fp):
    """
    The lookup table file is tab delimited and we're interested in the following columns:

    (0) (2)        (6)
    chr variant_id rsid_id_dbSNP147_GRCh37p13

    The variant_id is the unique GTEx identifier which is also found in the eQTL data
    files.

    arguments
        fp: filepath to the lookup table

    returns
        a dictionary mapping GTEx variant IDs to dbSNP identifiers (rsID)
    """

    table = []

    for df in pd.read_csv(fp, sep='\t', compression='infer', chunksize=4096):
        ## Rename ugly ass column name
        df = df.rename({'rs_id_dbSNP147_GRCh37p13': 'rsid'}, axis=1)
        ## Remove anything that doesn't have an rsID
        df = df[df.rsid != '.']

        table.extend(zip(df.variant_id, df.rsid))

    return dict([(v, r) for v, r in table])

def parse_eqtl_file(fp):
    """
    The eQTL file is tab delimited and we're interested in the following columns:

    (0)        (1)     (6)
    variant_id gene_id p-value

    arguments
        fp: eQTL filepath

    returns
        a data frame containing the parsed file. The frame only contains the 
        variant_id and gene_id columns.
    """

    df = pd.read_csv(fp, sep='\t', compression='infer')

    ## Remove the Ensembl gene version
    df['gene_id'] = df.gene_id.map(lambda g: g.split('.')[0])

    ## Remove anything w/ p < 0.05 (I don't think there are any but just in case)
    df = df[df.pval_nominal < 0.05]

    ## Rename p-value column
    df = df.rename(columns={'pval_nominal': 'p'})

    return df[['variant_id', 'gene_id', 'p']]

if __name__ == '__main__':
    from argparse import ArgumentParser

    ## cmd line shit
    usage = '%s [options] <snps> <merge> <output>' % argv[0]
    parse = ArgumentParser(usage=usage)

    parse.add_argument(
        'snps',
        nargs='?',
        help='tab delimited file containing SNP identifiers'
    )

    parse.add_argument(
        'merge',
        nargs='?',
        help='RsMergeArch table from NCBI'
    )

    parse.add_argument(
        'output',
        nargs='?',
        help='output file'
    )

    parse.add_argument(
        '-c',
        '--column',
        action='store',
        default='rsid',
        dest='column',
        type=str,
        help='column name in the <snps> file containing rsIDs (default = rsid)'
    )

    parse.add_argument(
        '--verbose',
        action='store_true',
        dest='verbose',
        help='Clutter your screen with output'
    )

    args = parse.parse_args()

    ## Add a console logger
    conlog = logging.StreamHandler()
    ## Set logging levels based on verbosity
    conlog.setLevel(logging.INFO if args.verbose else logging.ERROR)
    LOG.setLevel(logging.INFO if args.verbose else logging.ERROR)
    LOG.addHandler(conlog)

    if not args.snps:
        LOG.error('[!] You need to provide a file containing SNPs')
        LOG.error('')
        parse.print_help()
        exit(1)

    if not args.merge:
        LOG.error('[!] You need to provide a merge table from NCBI')
        LOG.error('')
        parse.print_help()
        exit(1)

    if not args.output:
        LOG.error('[!] You need to provide an output file')
        LOG.error('')
        parse.print_help()
        exit(1)

    LOG.info('[+] Reading merge table')

    ## Read in the merge table
    ## Table description can be found: 
    ## https://www.ncbi.nlm.nih.gov/projects/SNP/snp_db_table_description.cgi?t=RsMergeArch
    mtable = pd.read_csv(
        args.merge,
        sep='\t',
        compression='infer',
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
    
    # We really only need high and current for the mapping procedure
    mtable = mtable[['high', 'current']]

    LOG.info('[+] Reading SNPs dataset')

    ## Read in the snps dataset
    snps = pd.read_csv(args.snps, sep='\t')

    ## Strip the rs prefix, convert to int
    snps[args.column] = snps[args.column].str.strip(to_strip='rs')
    snps[args.column] = snps[args.column].astype(np.int64)

    LOG.info('[+] Updating old identifiers')

    ## Merge frames based on old SNP identifiers
    snps = snps.join(mtable.set_index('high'), on=args.column, how='left')

    criterion = snps.current.notna()

    ## Update IDs if necessary
    snps.loc[criterion, args.column] = snps.loc[criterion, 'current']

    snps = snps.drop(columns='current')

    ## Add the rs prefix back in
    snps[args.column] = 'rs' + snps[args.column].astype(np.int64).astype('str')

    snps.to_csv(
        args.output,
        sep='\t',
        index=False
    )

    LOG.info('[+] Done!')

