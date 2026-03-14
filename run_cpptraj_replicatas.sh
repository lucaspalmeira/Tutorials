#!/usr/bin/env bash

set -euo pipefail

# CONFIGURATION

# Topology file
PARM="step3_input.parm7"

# Last protein residue (REPLACE HERE)
# Example: if the protein goes from residue 1 to 412, use 412
PROT_LAST_RES="XXX"

# Protein mask
PROT_MASK=":1-${PROT_LAST_RES}"

# Ligand mask
LIG_MASK=":1CU,0CU"

# Number of clusters
NCLUSTERS=3

# To include solvent in H-bond analyses, set this to "yes"
INCLUDE_SOLVENT_HBOND="no"

# Solvent mask in the system
SOLVENT_MASK=":WAT"

# Replicates
REPLICATES=(1 2 3)

# =========================
# CHECKS
# =========================

if [[ ! -f "$PARM" ]]; then
    echo "ERROR: topology file not found: $PARM"
    exit 1
fi

if [[ "$PROT_LAST_RES" == "XXX" ]]; then
    echo "ERROR: replace PROT_LAST_RES=\"XXX\" with the last protein residue."
    exit 1
fi

if ! command -v cpptraj &>/dev/null; then
    echo "ERROR: cpptraj is not available in PATH."
    exit 1
fi

# REPLICATE LOOP

for i in "${REPLICATES[@]}"; do
    TRAJ="step5_production_${i}.nc"
    OUTDIR="replicate_analysis_${i}"

    echo "--------------------------------------"
    echo "Processing replicate ${i}"
    echo "Trajectory: ${TRAJ}"
    echo "Output: ${OUTDIR}"
    echo "--------------------------------------"

    if [[ ! -f "$TRAJ" ]]; then
        echo "WARNING: trajectory not found: $TRAJ"
        echo "Skipping replicate ${i}..."
        continue
    fi

    mkdir -p "$OUTDIR"

    # Trajectory centering
    cat > "${OUTDIR}/01_centering.in" << EOF
autoimage anchor ${PROT_MASK}
center ${PROT_MASK} mass origin
image origin center
rms first ${PROT_MASK}@CA
trajout ${OUTDIR}/step5_centered_${i}.nc netcdf
EOF

    cpptraj -p "$PARM" -y "$TRAJ" -i "${OUTDIR}/01_centering.in" \
        > "${OUTDIR}/01_centering.log" 2>&1

    # Ligand RMSD
    cat > "${OUTDIR}/02_rmsd_ligand.in" << EOF
rms Protein first ${PROT_MASK}@CA
rms Ligand first ${LIG_MASK} out ${OUTDIR}/rmsd_ligand.dat
EOF

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/02_rmsd_ligand.in" \
        > "${OUTDIR}/02_rmsd_ligand.log" 2>&1

    # Per-residue RMSF (Cα)
    cat > "${OUTDIR}/03_rmsf_ca.in" << EOF
rms first ${PROT_MASK}@CA
atomicfluct out ${OUTDIR}/rmsf_ca.dat ${PROT_MASK}@CA byres
EOF

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/03_rmsf_ca.in" \
        > "${OUTDIR}/03_rmsf_ca.log" 2>&1

    # Radius of gyration
    cat > "${OUTDIR}/04_radgyr.in" << EOF
radgyr ${PROT_MASK} out ${OUTDIR}/rg_protein.dat
EOF

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/04_radgyr.in" \
        > "${OUTDIR}/04_radgyr.log" 2>&1

    # Ligand–protein H-bonds
    if [[ "$INCLUDE_SOLVENT_HBOND" == "yes" ]]; then
        cat > "${OUTDIR}/05_hbond.in" << EOF
hbond HB out ${OUTDIR}/hbond_lig_prot.dat \
solventdonor ${SOLVENT_MASK} \
donormask ${LIG_MASK} \
acceptormask ${PROT_MASK}
EOF
    else
        cat > "${OUTDIR}/05_hbond.in" << EOF
hbond HB out ${OUTDIR}/hbond_lig_prot.dat \
donormask ${LIG_MASK} \
acceptormask ${PROT_MASK}
EOF
    fi

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/05_hbond.in" \
        > "${OUTDIR}/05_hbond.log" 2>&1

    # Ligand clustering
    cat > "${OUTDIR}/06_cluster.in" << EOF
rms first ${LIG_MASK}
cluster hieragglo clusters ${NCLUSTERS} linkage average \
    summary ${OUTDIR}/cluster_summary.dat \
    out ${OUTDIR}/cluster.dat \
    repout ${OUTDIR}/cluster_rep repfmt pdb
EOF

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/06_cluster.in" \
        > "${OUTDIR}/06_cluster.log" 2>&1

    echo "Replicate ${i} completed."
    echo
done

echo "All analyses completed."
