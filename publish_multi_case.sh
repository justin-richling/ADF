#!/usr/bin/env bash

set -euo pipefail

########################################
# Defaults
########################################
MODE=""
HOST=""
DEST=""
LIST_FILE=""
LOG_FILE="transfer.log"

########################################
# Helper Functions
########################################

usage() {
    cat <<EOF
Usage:
  $0 --mode local  --dest <directory> --list <file>
  $0 --mode remote --host <hostname> --dest <directory> --list <file>

Options:
  --mode     local or remote
  --host     remote hostname (required for remote mode)
  --dest     destination base directory
  --list     text file containing items (one per line)
  --log      optional log file (default: transfer.log)
EOF
    exit 1
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

########################################
# Argument Parsing
########################################

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --host) HOST="$2"; shift 2 ;;
        --dest) DEST="$2"; shift 2 ;;
        --list) LIST_FILE="$2"; shift 2 ;;
        --log)  LOG_FILE="$2"; shift 2 ;;
        *) usage ;;
    esac
done

########################################
# Validation
########################################

[[ -z "$MODE" ]] && usage
[[ -z "$DEST" ]] && usage
[[ -z "$LIST_FILE" ]] && usage
[[ ! -f "$LIST_FILE" ]] && { echo "List file not found."; exit 1; }

if [[ "$MODE" == "remote" ]]; then
    [[ -z "$HOST" ]] && usage
fi

########################################
# Choose Transfer Tool
########################################

if command_exists rsync; then
    TRANSFER_TOOL="rsync -avz"
else
    log "rsync not found, falling back to scp"
    TRANSFER_TOOL="scp -r"
fi

#TRANSFER_TOOL="scp -r"

########################################
# Main Loop
########################################

log "Starting transfer in $MODE mode"
log "Using tool: $TRANSFER_TOOL"

while IFS= read -r item; do

    [[ -z "$item" ]] && continue

    SOURCE="$item"

    if [[ "$MODE" == "local" ]]; then
        TARGET="${DEST}/"

        mkdir -p "$TARGET"

        log "Copying "$SOURCE/website/" "contents" -> $TARGET"
        $TRANSFER_TOOL "$SOURCE/website/." "$TARGET"

    else
        TARGET="${HOST}:${DEST}/"

        log "Copying "$SOURCE/website/" "contents" -> $TARGET"
        $TRANSFER_TOOL "$SOURCE/website/." "$TARGET"
    fi

done < "$LIST_FILE"

log "Transfer complete."
