#!/bin/bash
# Keepalive for watch_sn_a_ext_grid.py - restarts the watcher if it
# died before writing a DONE marker. Installed as a cron job.
BASE=/viper/ptmp1/sbrennan/SN2024pba
STATE=$BASE/logs/sn_a_ext_grid_watch.state
export PATH=/usr/bin:/bin:/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin

if pgrep -f "watch_sn_a_ext_grid.py" > /dev/null; then
    exit 0
fi
if grep -q "^DONE=" "$STATE" 2> /dev/null; then
    exit 0
fi
cd "$BASE" || exit 1
nohup python scripts/watch_sn_a_ext_grid.py \
    >> logs/sn_a_ext_grid_watch.nohup 2>&1 &
echo "$(date '+%F %T') keepalive: restarted watcher" \
    >> logs/sn_a_ext_grid_watch.log
