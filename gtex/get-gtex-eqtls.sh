#!/usr/bin/env bash

## file: get-gtex-eqtls.sh
## desc: Retrieves GTEx eQTL datasets.
## auth: TR

usage() {

    echo "usage: $0 [options]"
    echo ''
    echo 'Helper script to retrieve and pre-process GTEx eQTL datasets'
    echo ''
    echo 'Download options:'
    echo '  -f, --force  download eQTL data even if it already exists locally'
    echo ''
    echo 'Misc. options:'
    echo '  -h, --help   print this help message and exit'
    exit 0
}

## Load the config
[[ -f './config.sh' ]] && source './config.sh' || source '../config.sh'

## Cmd line options
force=""

while :; do
    case $1 in

        -f | --force)
            force=1
            ;;

        -h | -\? | --help)
            usage
            exit
            ;;

        --)
            shift; break;
            ;;

        -?*)
            log "WARN: unknown option (ignored): $1" >&2
            ;;

        *)
            break
    esac

    shift
done

## URL for the tissue-specific eQTL datasets
eqtl_url="https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz"
## URL for the eQTL lookup table used to map identifiers
table_url="https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz"
## URL for the eQTL annotations (used to derive tissue names/groups)
anno_url="https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt"

## Compressed output files
eqtl_gz="$DATA_DIR/gtex-eqtls-v7.tar.gz"
table_gz="$DATA_DIR/gtex-lookup-table.tsv.gz"
annos="$DATA_DIR/gtex-v7-annotations.tsv"

if [[ -n "$force" || ! -f "$eqtl_gz" ]]; then

    logger "Retrieving eQTL datasets"

    ## Options:
    ##  -q: quiet output
    ##  -O: save output to the specified file
    (wget -q "$eqtl_url" -O "$eqtl_gz") &
fi

if [[ -n "$force" || ! -f "$table_gz" ]]; then

    logger "Retrieving eQTL lookup table"

    (wget -q "$table_url" -O "$table_gz") &
fi

if [[ -n "$force" || ! -f "$annos" ]]; then

    logger "Retrieving GTEx annotations"

    (wget -q "$anno_url" -O "$annos") &
fi

wait

## Directory name in the tar'ed and zipped eQTL file
eqtl_dir="$DATA_DIR/GTEx_Analysis_v7_eQTL"

logger "Decompressing datasets"

tar -xf "$eqtl_gz" -C "$DATA_DIR"

## Rename the directory
mv "$eqtl_dir" "$DATA_DIR/gtex-eqtls"

