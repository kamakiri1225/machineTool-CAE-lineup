#!/bin/sh
cd "${0%/*}" || exit
chtMultiRegionFoam > log.stage1 2>&1 || exit 1
foamDictionary -entry PIMPLE.frozenFlow -add yes system/air/fvSolution > /dev/null
foamDictionary -entry endTime -set 21600 system/controlDict > /dev/null
foamDictionary -entry writeInterval -set 600 system/controlDict > /dev/null
foamDictionary -entry maxCo -set 1e5 system/controlDict > /dev/null
foamDictionary -entry maxDeltaT -set 30 system/controlDict > /dev/null
chtMultiRegionFoam > log.stage2 2>&1
touch post.foam
