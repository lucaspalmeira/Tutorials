#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Clusterização global entre 3 replicatas no cpptraj
# ============================================================
#
# Uso:
#   bash cluster_3_replicates.sh
#
# Opcional:
#   bash cluster_3_replicates.sh step3_input.parm7 combined_cluster_3reps
#
# Saída principal:
#   combined_cluster_3reps/most_representative_cluster_3reps.pdb
#
# ============================================================

PARM="${1:-step3_input.parm7}"
OUTDIR="${2:-combined_cluster_3reps}"

NCLUSTERS="${NCLUSTERS:-3}"

# Máscara usada para alinhar os frames antes da clusterização.
ALIGN_MASK="${ALIGN_MASK:-:1CU,0CU}"

# Máscara usada para calcular RMSD durante a clusterização.
# Por padrão usa apenas átomos pesados do ligante.
CLUSTER_MASK="${CLUSTER_MASK:-:1CU,0CU&!@H=}"

# Para incluir hidrogênios também na métrica de clusterização, rode:
#   CLUSTER_MASK=':1CU,0CU' bash cluster_3_replicates.sh
#
# Mas, para a maioria dos casos, é melhor clusterizar por átomos pesados
# e manter os hidrogênios apenas para análise posterior de H-bonds.

TRAJ1="replicate_analysis_1/step5_centered_1.nc"
TRAJ2="replicate_analysis_2/step5_centered_2.nc"
TRAJ3="replicate_analysis_3/step5_centered_3.nc"

mkdir -p "$OUTDIR"

echo "============================================================"
echo "Checando arquivos de entrada..."
echo "Topology:      $PARM"
echo "Traj 1:        $TRAJ1"
echo "Traj 2:        $TRAJ2"
echo "Traj 3:        $TRAJ3"
echo "Output:        $OUTDIR"
echo "N clusters:    $NCLUSTERS"
echo "Align mask:    $ALIGN_MASK"
echo "Cluster mask:  $CLUSTER_MASK"
echo "============================================================"

for f in "$PARM" "$TRAJ1" "$TRAJ2" "$TRAJ3"; do
    if [[ ! -f "$f" ]]; then
        echo "[ERRO] Arquivo não encontrado: $f"
        exit 1
    fi
done

CPPTRAJ_IN="$OUTDIR/combined_cluster_3reps.in"

cat > "$CPPTRAJ_IN" << EOF
parm $PARM

trajin $TRAJ1
trajin $TRAJ2
trajin $TRAJ3

# Alinhamento prévio.
rms first $ALIGN_MASK

# Clusterização global das três replicatas juntas.
cluster C0 \\
    hieragglo clusters $NCLUSTERS averagelinkage \\
    rms $CLUSTER_MASK \\
    out $OUTDIR/combined_cluster_vs_time.dat \\
    summary $OUTDIR/combined_cluster_summary.dat \\
    info $OUTDIR/combined_cluster_info.dat \\
    cpopvtime $OUTDIR/combined_cluster_population_vs_time.agr normframe \\
    repout $OUTDIR/combined_cluster_rep repfmt pdb \\
    avgout $OUTDIR/combined_cluster_avg avgfmt pdb

run
EOF

echo
echo "============================================================"
echo "Arquivo cpptraj gerado:"
echo "$CPPTRAJ_IN"
echo "============================================================"
cat "$CPPTRAJ_IN"

echo
echo "============================================================"
echo "Executando cpptraj..."
echo "============================================================"

cpptraj -i "$CPPTRAJ_IN" | tee "$OUTDIR/combined_cluster_3reps.log"

echo
echo "============================================================"
echo "Identificando o cluster mais populoso pelo summary..."
echo "============================================================"

python3 - << EOF
from pathlib import Path
import shutil

outdir = Path("$OUTDIR")
summary = outdir / "combined_cluster_summary.dat"

if not summary.exists():
    raise SystemExit(f"[ERRO] Arquivo não encontrado: {summary}")

clusters = []

with summary.open() as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("@"):
            continue

        parts = line.split()

        # Esperado:
        # Cluster Frames Frac AvgDist Stdev Centroid AvgCDist
        try:
            cid = int(parts[0])
            frames = int(parts[1])
            frac = float(parts[2])
        except Exception:
            continue

        clusters.append((cid, frames, frac))

if not clusters:
    raise SystemExit("[ERRO] Não consegui interpretar combined_cluster_summary.dat")

clusters_sorted = sorted(clusters, key=lambda x: x[1], reverse=True)
winner, nframes, frac = clusters_sorted[0]

print("População dos clusters:")
for cid, frames, f in sorted(clusters):
    print(f"  Cluster c{cid}: {frames} frames ({100*f:.2f}%)")

print()
print(f"Cluster mais representativo/global: c{winner}")
print(f"Número de frames: {nframes}")
print(f"Fração: {100*frac:.2f}%")

rep_pdb = outdir / f"combined_cluster_rep.c{winner}.pdb"

if not rep_pdb.exists():
    raise SystemExit(f"[ERRO] Não encontrei o PDB representativo esperado: {rep_pdb}")

final_pdb = outdir / "most_representative_cluster_3reps.pdb"
shutil.copy2(rep_pdb, final_pdb)

summary_out = outdir / "most_representative_cluster_3reps.txt"
with summary_out.open("w") as w:
    w.write("Global clustering across 3 replicates\\n")
    w.write("====================================\\n\\n")
    w.write(f"Most representative cluster: c{winner}\\n")
    w.write(f"Frames in cluster: {nframes}\\n")
    w.write(f"Fraction: {100*frac:.2f}%\\n")
    w.write(f"Representative PDB: {final_pdb.name}\\n\\n")
    w.write("Cluster populations:\\n")
    for cid, frames, f in sorted(clusters):
        w.write(f"c{cid}\\t{frames}\\t{100*f:.2f}%\\n")

print()
print(f"PDB final copiado para: {final_pdb}")
print(f"Resumo salvo em: {summary_out}")
EOF

echo
echo "============================================================"
echo "Finalizado!"
echo
echo "Arquivos principais:"
echo "  $OUTDIR/combined_cluster_summary.dat"
echo "  $OUTDIR/combined_cluster_info.dat"
echo "  $OUTDIR/combined_cluster_vs_time.dat"
echo "  $OUTDIR/most_representative_cluster_3reps.pdb"
echo "  $OUTDIR/most_representative_cluster_3reps.txt"
echo "============================================================"
