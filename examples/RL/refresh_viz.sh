#!/usr/bin/env bash
#
# refresh_viz.sh — regenerate the sheaf-RL documentation figures AND result numbers, then commit.
#
# By default it runs the visualizations ON THE CLUSTER (HiPerGator), where Julia 1.12 and the
# trained policies live: it syncs examples/RL up, runs every viz in one SLURM job, rsyncs the
# produced figures + metrics back into the repo, rewrites the result numbers in the docs from
# those metrics, and makes ONE local, subject-only commit. It NEVER pushes.
#
# Usage:
#   ./refresh_viz.sh              run on the cluster, bring figures+numbers back, commit
#   ./refresh_viz.sh --local      run the viz locally instead (fast; needs cache/*.jld2 present)
#   ./refresh_viz.sh --no-commit  do everything except the git commit
#   ./refresh_viz.sh --dry-run    print the plan only (no run, no copy, no commit)
#   ./refresh_viz.sh --only tv     restrict to figures whose key matches "tv"
#
# Numbers: each viz writes a `results_<tag>.toml` (analytic/oracle/rl). The docs mark every result
# with `<!--r:KEY-->value<!--/r-->` sentinels; this script rewrites the value between the sentinels
# from the collected metrics (see the KEY→metric map below).
set -euo pipefail

# --- locations ---------------------------------------------------------------
RL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$RL_DIR/../.." && pwd)"
CACHE="$RL_DIR/cache"
FIG="$REPO/docs/src/assets/figures/control"
DOCS_GLOB="$REPO/docs/src/control/sheaf_rl"          # where the number sentinels live

# cluster config (override via env)
HOST="${HPG_HOST:-hpg}"
REMOTE_REPO="${HPG_REPO:-/blue/fairbanksj/itaykadosh/CellularSheaves.jl}"
JULIA_REMOTE="${JULIA_REMOTE:-/apps/julia/1.12.0/bin/julia}"
JULIA_LOCAL="${JULIA:-julia}"
PROJECT_REL="examples/learned_mpc"                    # project that carries the viz deps (JLD2, Plots)

MODE=cluster; COMMIT=1; ONLY=""
while [ $# -gt 0 ]; do
  case "$1" in
    --local)     MODE=local ;;
    --dry-run)   MODE=dry ;;
    --no-commit) COMMIT=0 ;;
    --only)      ONLY="${2:-}"; shift ;;
    -h|--help)   sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "unknown flag: $1" >&2; exit 2 ;;
  esac
  shift
done

mkdir -p "$FIG"

# ---- the viz plan: key|script|env|produced-file(rel $CACHE)|doc figure(rel $FIG) ----
# '|' delimited because the env field itself contains '='. %C -> the cache dir at run time.
read -r -d '' PLAN <<'EOF' || true
tv|single_integrator.jl|ATV_SEEDS=9101 ATV_OUTDIR=%C/viz_tv|viz_tv/anim_tv_9101.gif|sheaf_tv_drift.gif
quad|quadrotor.jl|AQ_SEEDS=7001 AQ_OUTDIR=%C/viz_quad|viz_quad/anim_quad_7001.gif|sheaf_quad.gif
bc|bc_vs_analytic.jl|BCS_VIZ_SEEDS=9001 BCS_OUTDIR=%C/viz_bc|viz_bc/anim_bc_vs_analytic_9001.gif|sheaf_bc_clone.gif
multiagent_nom|multiagent.jl|AJ_TAG=nominal AJ_DRIFT=0 AJ_OUTDIR=%C/viz_multiagent|viz_multiagent/anim_nominal.gif|sheaf_multiagent_nominal.gif
multiagent|multiagent.jl|AJ_TAG=drift AJ_OUTDIR=%C/viz_multiagent|viz_multiagent/anim_drift.gif|sheaf_multiagent_drift.gif
multiagent_err|NONE||viz_multiagent/trkerr_drift.png|sheaf_multiagent_error.png
multiagent_3d|NONE||viz_multiagent/traj3d_drift.png|sheaf_multiagent_3d.png
multiagent_td|NONE||viz_multiagent/traj2d_drift.png|sheaf_multiagent_topdown.png
generalization|generalization.jl|EG_OUTDIR=%C/viz_multiagent|viz_multiagent/anim_generalization.gif|sheaf_multiagent_generalization.gif
hetero|heterogeneous.jl|AJ_TAG=drift AJ_OUTDIR=%C/viz_hetero|viz_hetero/anim_drift.gif|sheaf_hetero_drift.gif
hetero_err|NONE||viz_hetero/trkerr_drift.png|sheaf_hetero_error.png
hetero_3d|NONE||viz_hetero/anim3d_drift.gif|sheaf_hetero_3d.gif
ballcap|actuator_cap.jl|AJ_TAG=u15 AJ_OUTDIR=%C/viz_ballcap|viz_ballcap/anim_u15.gif|sheaf_ballcap_anim.gif
ballcap_err|NONE||viz_ballcap/trkerr_u15.png|sheaf_ballcap_error.png
EOF

# ---- KEY→(results file : metric) map for the number sentinels ----------------
# doc sentinel key = results_<tag>.toml metric.  d5=multiagent drift, d7=hetero, d8=ballcap,
# d2=tv, d3=quad (tv/quad emit lqr/analytic/oracle/rl as available).
read -r -d '' RESULTS_MAP <<'EOF' || true
d2_analytic=viz_tv/results_9101.toml:analytic
d2_oracle=viz_tv/results_9101.toml:oracle
d2_rl=viz_tv/results_9101.toml:rl
d3_lqr=viz_quad/results_7001.toml:analytic
d3_rl=viz_quad/results_7001.toml:rl
d5_analytic=viz_multiagent/results_drift.toml:analytic
d5_oracle=viz_multiagent/results_drift.toml:oracle
d5_rl=viz_multiagent/results_drift.toml:rl
d7_analytic=viz_hetero/results_drift.toml:analytic
d7_oracle=viz_hetero/results_drift.toml:oracle
d7_rl=viz_hetero/results_drift.toml:rl
d8_analytic=viz_ballcap/results_u15.toml:analytic
d8_rl=viz_ballcap/results_u15.toml:rl
EOF

want() { [ -z "$ONLY" ] && return 0; [[ "$1" == *"$ONLY"* ]]; }

# ---- 1. run the viz (cluster or local) --------------------------------------
emit_runner() {   # prints the shell commands to run all scripts, %C -> $1 (cache path)
  local cachep="$1"
  while IFS='|' read -r key script env src dst; do
    [ -z "$key" ] && continue
    want "$key" || continue
    [ "$script" = "NONE" ] && continue
    env="${env//%C/$cachep}"
    echo "echo '>> [$key] $script'; ( cd \"$VIZDIR\" && env $env \"\$JULIA\" --project=\"$PROJ\" \"$script\" ) || echo 'FAILED $key'"
  done <<< "$PLAN"
}

if [ "$MODE" = dry ]; then
  echo "# plan (dry run):"
  while IFS='|' read -r key script env src dst; do
    [ -z "$key" ] && continue; want "$key" || continue
    [ "$script" = NONE ] && echo "   place  $src -> $dst" || echo ">> [$key] $script  ($env)"
  done <<< "$PLAN"
  echo "# numbers rewritten from:"; while IFS= read -r m; do echo "   $m"; done <<< "$RESULTS_MAP"
  echo "(dry run — nothing executed)"; exit 0
fi

if [ "$MODE" = cluster ]; then
  echo "== cluster run on $HOST =="
  # sync scripts + lib + trained models up (into the cluster repo)
  rsync -az --delete -e "ssh -o BatchMode=yes" \
    --include='*/' --include='*.jl' --include='cache/***' --exclude='*' \
    "$RL_DIR/" "$HOST:$REMOTE_REPO/examples/RL/"
  VIZDIR="$REMOTE_REPO/examples/RL/viz"; PROJ="$REMOTE_REPO/$PROJECT_REL"
  RUNNER="$(JULIA='$JULIA' emit_runner "$REMOTE_REPO/examples/RL/cache")"
  # run everything in one SLURM job on a CPU partition
  ssh -o BatchMode=yes "$HOST" "cat > /tmp/sheaf_refresh_\$\$.sh" <<RSH
#!/bin/bash
#SBATCH --account=fairbanksj --qos=fairbanksj --partition=hpg-default
#SBATCH --cpus-per-task=4 --mem=16gb --time=00:40:00 --job-name=sheaf_viz
export GKSwstype=100
JULIA=$JULIA_REMOTE
$RUNNER
echo DONE_REFRESH
RSH
  JID=$(ssh -o BatchMode=yes "$HOST" "sbatch --parsable /tmp/sheaf_refresh_\$\$.sh" 2>/dev/null || true)
  echo "submitted SLURM job $JID; waiting..."
  ssh -o BatchMode=yes "$HOST" "while squeue -j $JID -h 2>/dev/null | grep -q .; do sleep 15; done"
  # bring the produced figures + metrics back
  rsync -az -e "ssh -o BatchMode=yes" "$HOST:$REMOTE_REPO/examples/RL/cache/" "$CACHE/"
elif [ "$MODE" = local ]; then
  echo "== local run =="
  VIZDIR="$RL_DIR/viz"; PROJ="$REPO/$PROJECT_REL"
  while IFS='|' read -r key script env src dst; do
    [ -z "$key" ] && continue; want "$key" || continue; [ "$script" = NONE ] && continue
    env="${env//%C/$CACHE}"
    echo ">> [$key] $script"
    ( cd "$VIZDIR" && env $env "$JULIA_LOCAL" --project="$PROJ" "$script" ) || echo "FAILED $key"
  done <<< "$PLAN"
fi

# ---- 2. place figures into the doc asset tree -------------------------------
while IFS='|' read -r key script env src dst; do
  [ -z "$key" ] && continue; want "$key" || continue
  if [ -f "$CACHE/$src" ]; then cp -f "$CACHE/$src" "$FIG/$dst"; echo "   placed $dst";
  else echo "   MISSING $CACHE/$src" >&2; fi
done <<< "$PLAN"

# ---- 3. rewrite the result numbers in the docs from the metrics -------------
toml_get() { grep -E "^$2 *=" "$CACHE/$1" 2>/dev/null | head -1 | sed -E 's/.*= *//' | tr -d ' '; }
if [ -d "$DOCS_GLOB" ]; then
  while IFS= read -r m; do
    [ -z "$m" ] && continue
    key="${m%%=*}"; filemetric="${m#*=}"; file="${filemetric%%:*}"; metric="${filemetric##*:}"
    val="$(toml_get "$file" "$metric")"
    [ -z "$val" ] && { echo "   (no metric for $key: $file:$metric)"; continue; }
    # rewrite <!--r:KEY-->OLD<!--/r--> -> <!--r:KEY-->val<!--/r-->  across the guide pages
    for md in "$DOCS_GLOB"/*.md; do
      [ -f "$md" ] || continue
      perl -0pi -e "s/(<!--r:${key}-->).*?(<!--\\/r-->)/\${1}${val}\${2}/g" "$md"
    done
    echo "   number $key = $val"
  done <<< "$RESULTS_MAP"
else
  echo "   (docs guide dir $DOCS_GLOB not found yet — skipping number rewrite)"
fi

# ---- 4. commit (local, subject-only, never push) ----------------------------
if [ "$COMMIT" -eq 1 ]; then
  cd "$REPO"; git add "$FIG" "$DOCS_GLOB" 2>/dev/null || git add "$FIG"
  if git diff --cached --quiet; then echo "no changes to commit";
  else git commit -m "Refresh sheaf-RL documentation figures and results"; echo "committed locally (not pushed)"; fi
fi
