#!/usr/bin/env python

## file: merge_snp_ids.py
## desc: Use the NCBI RsMergeArch table to merge old SNP identifiers into their
##       current versions.
## auth: TR

from __future__ import print_function
from sys import argv
import logging
import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)

if __name__ == '__main__':
    from argparse import ArgumentParser

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

