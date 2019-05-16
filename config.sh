#!/usr/bin/env bash

## file: config.sh
## desc: Contains configuration variables, such as data directories and commonly used
##       functions, for the shell scripts in this repo.

logger() { printf "[%s] %s\n" "$(date '+%Y.%m.%d %H:%M:%S')" "$*" >&2; }

## Retruns the directory that contains config.sh regardless of where the script is
## executed
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
## Source directory
SRC_DIR="$SELF_DIR/gtex"
## Data directory
DATA_DIR="$SELF_DIR/data"

DIRS=($DATA_DIR)

for d in "${DIRS[@]}"
do
    [[ ! -f "$d" ]] && mkdir -p "$d";
done

