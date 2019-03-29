#!/usr/bin/env bash

## file: get-merge-table.sh
## desc: Retrieves the reference SNP merge table from NCBI.
## auth: TR

usage() {

    echo "usage: $0 [options]"
    echo ''
    echo 'Helper script to retrieve the reference SNP merge table from NCBI.'
    echo ''
    echo 'Version options:'
    echo '  --v150       download dbSNP v. 150 data (default)'
    echo '  --v151       download dbSNP v. 151 data'
    echo '  --v152       download dbSNP v. 152 data'
    echo ''
    echo 'Download options:'
    echo '  -f, --force  download data even if it already exists locally'
    echo ''
    echo 'Misc. options:'
    echo '  -h, --help   print this help message and exit'
    exit 0
}

## Load the config
[[ -f './config.sh' ]] && source './config.sh' || source '../config.sh'

## Cmd line options
force=""
dbsnp=150

while :; do
    case $1 in

        --v150)
            ;;

        --v151)
            dbsnp=151
            ;;

        --v152)
            echo 'dbSNP v. 152 now reports merged identifiers in a JSON format.'
            echo 'Write your own parser for it.'
            exit 1
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

## URLs for the different dbSNP versions
URL150="ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/database/data/organism_data/RsMergeArch.bcp.gz"
URL151="ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/database/organism_data/RsMergeArch.bcp.gz"

if [[ "$dbsnp" -eq 150 ]]; then
    url="$URL150"
else
    url="$URL151"
fi

## Compressed output files
merge_gz="$DATA_DIR/dbsnp-merge-table.bcp.gz"

if [[ -n "$force" || ! -f "$merge_gz" ]]; then

    logger 'Retrieving NCBI merge table'

    ## Options:
    ##  -q: quiet output
    ##  -O: save output to the specified file
    (wget -q "$url" -O "$merge_gz") &

else
    logger 'Merge table exists, nothing to do.'
fi

