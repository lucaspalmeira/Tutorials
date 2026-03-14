#!/usr/bin/env bash

set -euo pipefail

# =========================
# CONFIGURAÇÕES
# =========================

# Topologia
PARM="step3_input.parm7"

# Último resíduo da proteína (SUBSTITUA AQUI)
# Exemplo: se a proteína vai do resíduo 1 ao 412, use 412
PROT_LAST_RES="XXX"

# Máscara da proteína
PROT_MASK=":1-${PROT_LAST_RES}"

# Máscara do ligante
LIG_MASK=":1CU,0CU"

# Número de clusters
NCLUSTERS=3

# Se quiser incluir solvente nas análises de H-bond, coloque "yes"
INCLUDE_SOLVENT_HBOND="no"

# Nome do solvente no sistema
SOLVENT_MASK=":WAT"

# Replicatas
REPLICAS=(1 2 3)

# CHECAGENS

if [[ ! -f "$PARM" ]]; then
    echo "ERRO: arquivo de topologia não encontrado: $PARM"
    exit 1
fi

if [[ "$PROT_LAST_RES" == "XXX" ]]; then
    echo "ERRO: substitua PROT_LAST_RES=\"XXX\" pelo último resíduo da proteína."
    exit 1
fi

if ! command -v cpptraj &>/dev/null; then
    echo "ERRO: cpptraj não está no PATH."
    exit 1
fi

# LOOP DAS REPLICATAS

for i in "${REPLICAS[@]}"; do
    TRAJ="step5_production_${i}.nc"
    OUTDIR="analise_replica_${i}"

    echo "======================================"
    echo "Processando replicata ${i}"
    echo "Trajetória: ${TRAJ}"
    echo "Saída: ${OUTDIR}"
    echo "======================================"

    if [[ ! -f "$TRAJ" ]]; then
        echo "AVISO: trajetória não encontrada: $TRAJ"
        echo "Pulando replicata ${i}..."
        continue
    fi

    mkdir -p "$OUTDIR"

    # Centralização da trajetória
    
    cat > "${OUTDIR}/01_centering.in" << EOF
autoimage anchor ${PROT_MASK}
center ${PROT_MASK} mass origin
image origin center
rms first ${PROT_MASK}@CA
trajout ${OUTDIR}/step5_centered_${i}.nc netcdf
EOF

    cpptraj -p "$PARM" -y "$TRAJ" -i "${OUTDIR}/01_centering.in" \
        > "${OUTDIR}/01_centering.log" 2>&1

    # RMSD do ligante
    
    cat > "${OUTDIR}/02_rmsd_ligand.in" << EOF
rms Protein first ${PROT_MASK}@CA
rms Ligand first ${LIG_MASK} out ${OUTDIR}/rmsd_ligand.dat
EOF

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/02_rmsd_ligand.in" \
        > "${OUTDIR}/02_rmsd_ligand.log" 2>&1

    
    # RMSF por resíduo (Cα)
    
    cat > "${OUTDIR}/03_rmsf_ca.in" << EOF
rms first ${PROT_MASK}@CA
atomicfluct out ${OUTDIR}/rmsf_ca.dat ${PROT_MASK}@CA byres
EOF

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/03_rmsf_ca.in" \
        > "${OUTDIR}/03_rmsf_ca.log" 2>&1

    # Raio de giro
    
    cat > "${OUTDIR}/04_radgyr.in" << EOF
radgyr ${PROT_MASK} out ${OUTDIR}/rg_protein.dat
EOF

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/04_radgyr.in" \
        > "${OUTDIR}/04_radgyr.log" 2>&1

    # H-bonds ligante–proteína
    
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

    # Clustering do ligante
    
    cat > "${OUTDIR}/06_cluster.in" << EOF
rms first ${LIG_MASK}
cluster hieragglo clusters ${NCLUSTERS} linkage average \
    summary ${OUTDIR}/cluster_summary.dat \
    out ${OUTDIR}/cluster.dat \
    repout ${OUTDIR}/cluster_rep repfmt pdb
EOF

    cpptraj -p "$PARM" -y "${OUTDIR}/step5_centered_${i}.nc" -i "${OUTDIR}/06_cluster.in" \
        > "${OUTDIR}/06_cluster.log" 2>&1

    echo "Replicata ${i} finalizada."
    echo
done

echo "Todas as análises concluídas."
