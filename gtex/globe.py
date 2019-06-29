#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## file: globe.py
## desc: GTEx eQTL pipeline global variables: directories, filepaths, and dataset URLs.
## auth: TR

import logging
import os

_logger = logging.getLogger(__name__)


## URLs ##

## URL for the tissue-specific eQTL datasets
_url_eqtls = 'https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz'

## URL for the eQTL lookup table used to map identifiers
_url_table = 'https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz'

## URL for the eQTL annotations (used to derive tissue names/groups)
_url_annotations = 'https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt'

## URLs for the NCBI dbSNP merge table. NCBI now reports merged SNPs in JSON format for
## all dbSNP versions >= 152. So write your own parser for that.
_url_dbsnp150 = "http://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/database/data/organism_data/RsMergeArch.bcp.gz"
_url_dbsnp151 = "http://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/database/organism_data/RsMergeArch.bcp.gz"

## UCSC tools and data
_url_ucsc_liftover = 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/liftOver'
_url_ucsc_hg38_chain = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'

## Output directories ##
_dir_data = 'data'

## Tools
_dir_data_tools = 'data/tools'

## Raw datasets
_dir_data_raw = 'data/raw'

## Complete, processed datasets
_dir_data_processed = 'data/processed'

## Complete, processed, and lifted datasets
_dir_data_lifted = 'data/lifted'

## Extracted eQTLs
_dir_eqtls = os.path.join(_dir_data_raw, 'GTEx_Analysis_v7_eQTL')


## Output files ##

## GTEx eQTLs
_fp_compressed_eqtls = os.path.join(_dir_data_raw, 'gtex-eqtls.tar.gz')

## GTEx lookup table
_fp_compressed_lookup_table = os.path.join(_dir_data_raw, 'gtex-lookup-table.tsv.gz')
_fp_lookup_table = os.path.join(_dir_data_raw, 'gtex-lookup-table.tsv')

## GTEx tissue/eQTL annotations
_fp_annotations = os.path.join(_dir_data_raw, 'gtex-annotations.tsv')

## NCBI dbSNP merge table
_fp_compressed_dbsnp_table = os.path.join(_dir_data_raw, 'dbsnp-merge-table.tsv.gz')
_fp_dbsnp_table = os.path.join(_dir_data_raw, 'dbsnp-merge-table.tsv')

## eQTL stats
_fp_eqtl_stats = os.path.join(_dir_data, 'gtex-eqtl-stats.tsv')

## liftOver executable
_exe_liftover = os.path.join(_dir_data_tools, 'liftOver')

## Zipped genome chains
_fp_hg38_chain_gz = os.path.join(_dir_data_tools, 'hg19-hg38.chain.gz')

## Unzipped genome chains
_fp_hg38_chain = os.path.join(_dir_data_tools, 'hg19-hg38.chain')


## In case these don't exist
try:
    os.makedirs(_dir_data_tools, exist_ok=True)
    os.makedirs(_dir_data_raw, exist_ok=True)
    os.makedirs(_dir_data_processed, exist_ok=True)
    os.makedirs(_dir_data_lifted, exist_ok=True)

except OSError as e:
    _logger.error('Failed to create data directories: %s', e)
    exit(1)

