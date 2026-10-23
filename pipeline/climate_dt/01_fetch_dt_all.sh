#!/bin/sh
# Run all Climate DT point extractions, 6 at a time (DestinE upgraded access, ~/.polytopeapirc).
# Usage (from the repository root, inside the pixi environment): sh pipeline/climate_dt/01_fetch_dt_all.sh work/climate_dt
ROOT=$1
{ for m in IFS-NEMO IFS-FESOM; do
    for c in 3490355 3490408 3490451 3490362 3490357; do echo "$m avg_tprate $c"; done
    for c in 3490355 3490451 3490362; do echo "$m 2t $c"; done
    echo "$m avg_ie 3490355"
  done; } | xargs -P 6 -L 1 sh -c "python -I pipeline/climate_dt/01_fetch_dt.py \"$ROOT\" \$0 \$1 \$2"
