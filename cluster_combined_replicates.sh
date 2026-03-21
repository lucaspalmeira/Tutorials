#!/usr/bin/env bash

set -euo pipefail

# =========================================================
# COMBINED CLUSTERING OF 3 MD REPLICATES WITH CPPTRAJ
# =========================================================

# =========================
# HELPER: LOGGING
# =========================

log() {
    echo "$@" >&2
}

die() {
    echo "ERROR: $@" >&2
    exit 1
}

# =========================
# USAGE
# =========================

usage() {
    cat >&2 <<EOF
Usage:
  bash $0 <active_site|contacts>

Examples:
  bash $0 active_site
  bash $0 contacts
EOF
    exit 1
}

# =========================
# CLUSTERING MODE ARGUMENT
# =========================

MODE="${1:-}"

[[ -n "$MODE" ]] || usage
[[ "$MODE" == "active_site" || "$MODE" == "contacts" ]] || \
    die "MODE must be either 'active_site' or 'contacts'"

# =========================
# USER SETTINGS
# =========================

# Topology file
PARM="step3_input.parm7"

# Production trajectories
TRAJ1="replicate_analysis_1/step5_centered_1.nc"
TRAJ2="replicate_analysis_2/step5_centered_2.nc"
TRAJ3="replicate_analysis_3/step5_centered_3.nc"

# Frame reading setup for trajin
TRAJ_START=1
TRAJ_STOP=last
TRAJ_STRIDE=10

# Ligand residue names/numbers in AMBER mask syntax
LIG_MASK=":6CU,0CU"

# Active-site residues
ACTIVE_SITE_RES=":13,146,205"

# Contact residues
CONTACT_RES=":13,146,205"

# Solvent/ions to remove before clustering
STRIP_MASK=":WAT,Na+,Cl-"

# Clustering method settings
CLUSTER_EPSILON="2.0"
NCLUSTERS="10"
SIEVE="10"

# Output directory
OUTDIR="combined_clustering"

# Prefix used in output file names
PREFIX="combined"

# =========================
# CHECKS
# =========================

command -v cpptraj >/dev/null 2>&1 || die "cpptraj was not found in PATH"

[[ -f "$PARM" ]] || die "topology file not found: $PARM"

for f in "$TRAJ1" "$TRAJ2" "$TRAJ3"; do
    [[ -f "$f" ]] || die "trajectory file not found: $f"
done

mkdir -p "$OUTDIR"

# =========================
# HELPER FUNCTIONS
# =========================

count_frames() {
    local traj="$1"
    local tmp_input
    local cpptraj_output
    local nframes

    tmp_input=$(mktemp)

    cat > "$tmp_input" <<EOF
parm $PARM
trajin $traj $TRAJ_START $TRAJ_STOP $TRAJ_STRIDE
run
quit
EOF

    if ! cpptraj_output=$(cpptraj -i "$tmp_input" 2>&1); then
        rm -f "$tmp_input"
        log "ERROR while reading trajectory: $traj"
        log "---- cpptraj output begin ----"
        log "$cpptraj_output"
        log "---- cpptraj output end ----"
        return 1
    fi

    rm -f "$tmp_input"

    nframes=$(echo "$cpptraj_output" | grep -Eo 'Frames:[[:space:]]*[0-9]+' | head -n1 | grep -Eo '[0-9]+' || true)

    if [[ -z "${nframes:-}" ]]; then
        nframes=$(echo "$cpptraj_output" | grep -Eo 'Read[[:space:]]+[0-9]+[[:space:]]+frames' | head -n1 | grep -Eo '[0-9]+' | head -n1 || true)
    fi

    if [[ -z "${nframes:-}" ]]; then
        nframes=$(echo "$cpptraj_output" | awk '
            /Frames:/ {
                for (i = 1; i <= NF; i++) {
                    if ($i == "Frames:") {
                        print $(i+1)
                        exit
                    }
                }
            }
        ' || true)
    fi

    if [[ -z "${nframes:-}" ]]; then
        log "ERROR: could not determine number of frames for $traj"
        log "---- cpptraj output begin ----"
        log "$cpptraj_output"
        log "---- cpptraj output end ----"
        return 1
    fi

    echo "$nframes"
}

validate_mask() {
    local label="$1"
    local mask="$2"
    local tmp_input
    local cpptraj_output

    tmp_input=$(mktemp)

    cat > "$tmp_input" <<EOF
parm $PARM
select "$mask"
run
quit
EOF

    if ! cpptraj_output=$(cpptraj -i "$tmp_input" 2>&1); then
        rm -f "$tmp_input"
        log "ERROR while validating $label mask: $mask"
        log "---- cpptraj output begin ----"
        log "$cpptraj_output"
        log "---- cpptraj output end ----"
        return 1
    fi

    rm -f "$tmp_input"

    if ! echo "$cpptraj_output" | grep -qi "Selected="; then
        log "ERROR: could not validate $label mask: $mask"
        log "---- cpptraj output begin ----"
        log "$cpptraj_output"
        log "---- cpptraj output end ----"
        return 1
    fi

    log "Mask OK [$label]: $mask"
}

build_cluster_mask() {
    if [[ "$MODE" == "active_site" ]]; then
        # Ligand + active-site residues, heavy atoms only
        echo "(((${LIG_MASK})|(${ACTIVE_SITE_RES}))&!@H=)"
    else
        # Contacts mode:
        # 1) ligand heavy atoms
        # 2) CA atoms of contact residues
        # 3) side-chain heavy atoms of contact residues
        echo "(((${LIG_MASK})&!@H=)|(${CONTACT_RES}@CA)|(((${CONTACT_RES})&!@H=)&!@N,CA,C,O,OXT))"
    fi
}

# =========================
# VALIDATE MASKS
# =========================

log "Validating masks..."
validate_mask "ligand" "$LIG_MASK"
validate_mask "active-site residues" "$ACTIVE_SITE_RES"
validate_mask "contact residues" "$CONTACT_RES"

# =========================
# FRAME COUNT FOR SPLITFRAME
# =========================

log ""
log "Counting frames after trajin selection..."

N1=$(count_frames "$TRAJ1") || die "failed to count frames for $TRAJ1"
N2=$(count_frames "$TRAJ2") || die "failed to count frames for $TRAJ2"
N3=$(count_frames "$TRAJ3") || die "failed to count frames for $TRAJ3"

SPLIT1="$N1"
SPLIT2=$((N1 + N2))
TOTAL=$((N1 + N2 + N3))

log "Frames:"
log "  Replica 1: $N1"
log "  Replica 2: $N2"
log "  Replica 3: $N3"
log "  Total:     $TOTAL"
log "splitframe = ${SPLIT1},${SPLIT2}"

# =========================
# BUILD MASK
# =========================

CLUSTER_MASK=$(build_cluster_mask)

log ""
log "Clustering mode: $MODE"
log "Cluster mask: $CLUSTER_MASK"

# =========================
# OPTIONAL: VALIDATE FINAL CLUSTER MASK
# =========================

log ""
log "Validating final cluster mask..."
validate_mask "cluster mask" "$CLUSTER_MASK"

# =========================
# PREPARE CPPTRAJ INPUT
# =========================

CPPTRAJ_INPUT="${OUTDIR}/cluster_${MODE}.in"
CPPTRAJ_LOG="${OUTDIR}/cluster_${MODE}.log"
OUTPREFIX="${OUTDIR}/${PREFIX}_${MODE}"

cat > "$CPPTRAJ_INPUT" <<EOF
parm $PARM

trajin $TRAJ1 $TRAJ_START $TRAJ_STOP $TRAJ_STRIDE
trajin $TRAJ2 $TRAJ_START $TRAJ_STOP $TRAJ_STRIDE
trajin $TRAJ3 $TRAJ_START $TRAJ_STOP $TRAJ_STRIDE

autoimage
strip $STRIP_MASK

cluster C0 \
    hieragglo epsilon $CLUSTER_EPSILON clusters $NCLUSTERS averagelinkage \
    rms "$CLUSTER_MASK" \
    sieve $SIEVE random \
    out ${OUTPREFIX}.cnumvtime.dat \
    summary ${OUTPREFIX}.summary.dat \
    summarysplit ${OUTPREFIX}.split.dat splitframe ${SPLIT1},${SPLIT2} \
    info ${OUTPREFIX}.info.dat \
    cpopvtime ${OUTPREFIX}.cpopvtime.agr normframe \
    repout ${OUTPREFIX}.rep repfmt pdb \
    avgout ${OUTPREFIX}.avg avgfmt pdb

run
quit
EOF

# =========================
# RUN CPPTRAJ
# =========================

log ""
log "Running cpptraj..."
cpptraj -i "$CPPTRAJ_INPUT" | tee "$CPPTRAJ_LOG"

# =========================
# FINAL REPRESENTATIVE PDB
# =========================

REP_PDB="${OUTPREFIX}.rep.c0.pdb"
FINAL_PDB="${OUTDIR}/representative_global_${MODE}.pdb"

[[ -f "$REP_PDB" ]] || die "representative PDB for cluster 0 was not generated: $REP_PDB"

cp "$REP_PDB" "$FINAL_PDB"

# =========================
# DONE
# =========================

log ""
log "Done."
log "Main outputs:"
log "  cpptraj input:              $CPPTRAJ_INPUT"
log "  cpptraj log:                $CPPTRAJ_LOG"
log "  cluster summary:            ${OUTPREFIX}.summary.dat"
log "  cluster split summary:      ${OUTPREFIX}.split.dat"
log "  cluster info:               ${OUTPREFIX}.info.dat"
log "  representative PDB (c0):    $REP_PDB"
log "  final copied PDB:           $FINAL_PDB"
