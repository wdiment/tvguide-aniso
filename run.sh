#!/bin/bash

### Runner script for TVGuide-Aniso
### Script will check for a local tvguide environment and the mamba/conda commands
### If not found, will initiate build script
### If found, will start tvguide-aniso with custom gunicorn settings
### Gunicorn settings can be modified in order to fit system settings as necessary


### This script lives in the repo root - always run from there, whatever the folder is named
cd "$(dirname "${BASH_SOURCE[0]}")"
export MINIFORGE_LOCATION=$PWD


bash $MINIFORGE_LOCATION/build.sh
if [[ -z `command -v conda` || -z `command -v mamba`  ]]; then
    MINIFORGE_PREFIX="$MINIFORGE_LOCATION/.miniforge3"
    if [ -d $MINIFORGE_PREFIX ]; then
        source "$MINIFORGE_PREFIX/etc/profile.d/conda.sh"
    else
        echo "Mamba/Conda commands not found - exiting"
        exit 1
    fi
fi

### Optional local settings (gitignored): MAPBOX_TOKEN, TVGUIDE_BIND, etc.
if [ -f "$MINIFORGE_LOCATION/.env" ]; then
    set -a
    source "$MINIFORGE_LOCATION/.env"
    set +a
fi

echo "Starting TVGuide-Aniso "
eval "$(conda shell.bash hook)"
conda activate $MINIFORGE_LOCATION/tvguide

### Stop a previous instance of this app only - not every gunicorn on the machine
PIDFILE="$MINIFORGE_LOCATION/tvguide.pid"
if [ -f "$PIDFILE" ]; then
    kill "$(cat "$PIDFILE")" 2>/dev/null
    rm -f "$PIDFILE"
fi
### Kill any other previous gunicorn instances
pkill -f "gunicorn.*tvguide:server" 2>/dev/null
sleep 1

gunicorn -w 3 -t 6 --timeout 90 -b "${TVGUIDE_BIND:-0.0.0.0:8000}" --pid "$PIDFILE" tvguide:server
