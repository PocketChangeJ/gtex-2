#!/usr/bin/env bash

## file: process-brain-eqtls.sh
## desc: Sample pipeline to retrieve and process brain eQTLs from GTEx.
## auth: TR

## Load the config
[[ -f './config.sh' ]] && source './config.sh' || source '../config.sh'

usage() {

    echo "usage: $0 [options]"
    echo ''
    echo 'Retrieve, process, and format brain eQTLs from GTEx.'
    echo ''
    echo 'Misc. options:'
    echo '  -h, --help   print this help message and exit'
    exit 0
}

while :; do
    case $1 in

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

## Filpaths to the eQTL directory, lookup table, annotations, merge table,
## unmapped eQTLs, and processed output
eqtls="$DATA_DIR/gtex-eqtls"
lookup="$DATA_DIR/gtex-lookup-table.tsv.gz" 
annos="$DATA_DIR/gtex-v7-annotations.tsv"
merge="$DATA_DIR/dbsnp-merge-table.bcp.gz"
nomap="$DATA_DIR/unmapped-gtex-brain-eqtls.tsv"
output="$DATA_DIR/gtex-brain-eqtls.tsv"

## Retrieve GTEx eQTLs
"$SRC_DIR/get-gtex-eqtls.sh"

## Process all the brain eQTLs (13 tissues in total)
python "$SRC_DIR/process_gtex_eqtls.py" --verbose -d -f brain -a "$annos" -u "$nomap" "$eqtls" "$lookup" "$output"

## Retrieve the NCBI merge table to convert old reference SNP identifiers
"$SRC_DIR/get-merge-table.sh"

## Update old IDs
python "$SRC_DIR/merge_snp_ids.py" --verbose "$output" "$merge" "${output%.*}-dbsnp150.tsv"

