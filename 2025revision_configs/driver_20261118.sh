#!/bin/bash
# Sequential restart (user 2026-11-18): encore_edgecases_K562 -> encore_edgecases_HepG2; stop if <5T free on ps-yeolab3
B=/tscc/nfs/home/bay001/projects/charlene_paper_20250521/permanent_data/skipper2025revision
for r in encore_edgecases_K562 encore_edgecases_HepG2; do
    avail=$(df -B1 --output=avail /tscc/projects/ps-yeolab3 | tail -1)
    echo "$(date) df: $(df -h /tscc/projects/ps-yeolab3 | tail -1)"
    [ $avail -ge 5000000000000 ] || { echo "$(date) STOP: <5T free on ps-yeolab3"; exit 1; }
    echo "$(date) START $r"
    $B/$r/run_skipper.sh > $B/$r/supervise.log 2>&1; rc=$?
    echo "$(date) END $r rc=$rc"
done
echo "$(date) DRIVER DONE"
