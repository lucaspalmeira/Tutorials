#!/usr/bin/env bash

set -e

# Minim
pmemd.cuda -O \
-i step4.0_minimization.mdin \
-p step3_input.parm7 \
-c step3_input.rst7 \
-o step4.0_minimization.mdout \
-r step4.0_minimization.rst7 \
-inf step4.0_minimization.mdinfo \
-ref step3_input.rst7

# Equil
pmemd.cuda -O \
-i step4.1_equilibration.mdin \
-p step3_input.parm7 \
-c step4.0_minimization.rst7 \
-o step4.1_equilibration.mdout \
-r step4.1_equilibration.rst7 \
-inf step4.1_equilibration.mdinfo \
-ref step3_input.rst7 \
-x step4.1_equilibration.nc


# Prod (replicates)
for i in 1 2 3
do
    echo "Run replicate $i"

    pmemd.cuda -O \
    -i step5_production.mdin \
    -p step3_input.parm7 \
    -c step4.1_equilibration.rst7 \
    -o step5_production_${i}.mdout \
    -r step5_production_${i}.rst7 \
    -inf step5_production_${i}.mdinfo \
    -x step5_production_${i}.nc

done
