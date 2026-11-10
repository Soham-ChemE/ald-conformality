#!/bin/sh
# Run a JSL script in JMP Student Edition 19 from the shell and wait for it to finish.
# JMP opens .jsl files in the editor unless the first line is "//!", so a runnable copy is made with that marker.
# The script itself must write its outputs to files and end with Quit( No Save ) (or leave JMP open if you want).
# Usage: sh scripts/run_jsl.sh path/to/script.jsl [seconds to wait, default 120]
set -e
SRC="$1"; WAIT="${2:-120}"
APP="/Applications/JMP Student Edition 19.app"
TMP="$(mktemp -t jsl).jsl"
{ echo '//!'; cat "$SRC"; } > "$TMP"
open -g -W -a "$APP" "$TMP" &
PID=$!
i=0; while kill -0 $PID 2>/dev/null && [ $i -lt "$WAIT" ]; do sleep 1; i=$((i+1)); done
rm -f "$TMP"
