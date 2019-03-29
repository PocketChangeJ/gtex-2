#!/usr/bin/env python

## file: process_gtex_eqtls.py
## desc: Parses GTEx eQTL data sets and formats them for later use with GeneWeaver.
## auth: TR

from __future__ import print_function
from sys import argv
import logging
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
        a data frame containing the parsed file. The frame only contains the variant_id
        and gene_id columns.
    """

    df = pd.read_csv(fp, sep='\t', compression='infer')

    ## Remove the Ensembl gene version
    df.gene_id = df.gene_id.map(lambda g: g.split('.')[0])

    ## Remove anything w/ p < 0.05 (I don't think there are any but just in case)
    df = df[df.pval_nominal < 0.05]

    ## Rename p-value column
    df = df.rename(columns={'pval_nominal': 'p'})

    return df[['variant_id', 'gene_id', 'p']]

if __name__ == '__main__':
    from argparse import ArgumentParser

    ## cmd line shit
    usage = '%s [options] <eqtls> <table> <output>' % argv[0]
    parse = ArgumentParser(usage=usage)

    parse.add_argument(
        'eqtls',
        nargs='?',
        help='file containing significant GTEx eQTL variant-gene pairs'
    )

    parse.add_argument(
        'table',
        nargs='?',
        help='GTEx eQTL lookup table'
    )

    parse.add_argument(
        'output',
        nargs='?',
        help='output file'
    )

    parse.add_argument(
        '-a',
        '--annotations',
        action='store',
        dest='annotations',
        default='',
        metavar='<file>',
        type=str,
        help='use the given GTEx annotations to annotate tissues'
    )

    parse.add_argument(
        '-d',
        '--dir',
        action='store_true',
        dest='dir',
        help='use a directory of eQTL datasets instead of a single file'
    )

    parse.add_argument(
        '-f',
        '--filter',
        action='store',
        dest='filter',
        metavar='<list>',
        help='filter eQTLs based on tissue groups (requires -a/--annotations)'
    )

    parse.add_argument(
        '-t',
        '--tissue',
        action='store',
        dest='tissue',
        default='unknown',
        metavar='<tissue>',
        help='annotate eQTLs using the given tissue'
    )

    parse.add_argument(
        '-u',
        '--unmapped',
        action='store',
        dest='unmapped',
        default='unmapped-gtex-eqtls.tsv',
        metavar='<file>',
        help='store unmapped eQTLs at the given filepath'
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

    if not args.eqtls:
        LOG.error('[!] You need to provide an eQTL dataset')
        LOG.error('')
        parse.print_help()
        exit(1)

    if not args.table:
        LOG.error('[!] You need to provide a GTEx lookup table')
        LOG.error('')
        parse.print_help()
        exit(1)

    if not args.output:
        LOG.error('[!] You need to provide an output file')
        LOG.error('')
        parse.print_help()
        exit(1)

    ## Get eQTL files
    if args.dir:
        LOG.info('[+] Enemurating eQTL datasets')

        eqtl_files = find_eqtl_files(args.eqtls)
    else:
        eqtl_files = [args.eqtls]

    if not eqtl_files:
        LOG.error('[!] No eQTL datasets were found')
        exit(1)

    LOG.info('[+] Parsing lookup table')

    lookup = parse_lookup_file(args.table)

    ## Parse annotations if necessary
    if args.annotations:
        LOG.info('[+] Parsing annotations')

        annotations = parse_annotations_file(args.annotations)
    else:
        annotations = pd.DataFrame()

    eqtl_list = []
    unmapped_list = []

    for fl in eqtl_files:
        LOG.info('[+] Parsing %s', fl)

        eqtls = parse_eqtl_file(fl)

        eqtls['tissue'] = args.tissue
        eqtls['tissue_group'] = args.tissue

        ## Add a tissue if possible
        if not annotations.empty:

            for row in annotations.itertuples():

                ## Try to match to an annotation
                if re.search(row.tissue_str, fl, flags=re.IGNORECASE):
                    eqtls['tissue'] = row.tissue_name
                    eqtls['tissue_group'] = row.tissue_group

                    break
            else:
                LOG.warn('[-] Could not discern tissue type for %s', fl)

            ## Filter if necessary
            if args.filter:
                eqtls = eqtls[eqtls.tissue_group.str.contains(args.filter, case=False)]

        ## Map GTEx eQTL IDs to dbSNP references
        eqtls['rsid'] = eqtls.variant_id.map(lambda v: lookup.get(v))

        ## Seprate and remove eQTLs that didn't map to an rsID
        unmapped = eqtls[eqtls.rsid.isnull()]
        eqtls = eqtls[eqtls.rsid.notnull()]

        eqtl_list.append(eqtls)
        unmapped_list.append(unmapped)

    ## Merge dataframes
    eqtl_list = pd.concat(eqtl_list)
    unmapped_list = pd.concat(unmapped_list)

    LOG.info('[+] Saving eQTL data')

    eqtl_list.to_csv(
        args.output,
        sep='\t',
        columns=['variant_id', 'rsid', 'gene_id', 'p', 'tissue', 'tissue_group'],
        index=False
    )

    ## If there are unmapped eQTLs, save those as well
    if not unmapped_list.empty:
        LOG.warn('[!] Some eQTLs were not mapped to dbSNP reference identifiers')

        unmapped_list.to_csv(
            args.unmapped,
            sep='\t',
            columns=['variant_id', 'gene_id', 'p', 'tissue', 'tissue_group'],
            index=False
        )

    LOG.info('[+] Done!')

