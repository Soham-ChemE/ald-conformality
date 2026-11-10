#!/bin/sh
# Run inside the OpenFOAM env: for each temperature case run steady then pulse (direct calls, no bash subshell).
cd "$(dirname "$0")/temperature"
for C in T200_blanket T200_patterned T300_blanket T300_patterned T400_blanket T400_patterned; do
  ( cd $C/steady && blockMesh > log.blockMesh 2>&1 && simpleFoam > log.simpleFoam 2>&1 )
  LAST=$(ls -d $C/steady/[1-9]* | sort -t/ -k3 -g | tail -1)
  rm -rf $C/pulse/constant/polyMesh; cp -r $C/steady/constant/polyMesh $C/pulse/constant/; cp $LAST/U $LAST/phi $C/pulse/0/
  ( cd $C/pulse && scalarTransportFoam > log.scalarTransportFoam 2>&1; echo "$C exit $? steps $(grep -c '^Time = ' log.scalarTransportFoam)" )
done
