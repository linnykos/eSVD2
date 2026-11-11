#!/usr/bin/env bash
# Log one experiment of an autoresearch loop and move the branch accordingly.
#
# Install: copy this file to <loop>/harness/record_result.sh and write <loop>/harness/loop.conf beside it with
#   TRIAL_FILE=<path of the file on trial, relative to the repository root>
#   DIRECTION=higher|lower        the score's direction
#   EXTRA_FIELDS="<block keys>"   block keys copied into results.tsv as columns, space-separated; may be empty
#   SCORER="<command>"            how the scorer is run, for messages only (e.g. "Rscript autoresearch_x/harness/evaluate.R" or
#                                 "python autoresearch_x/harness/evaluate.py")
# The loop folder is the parent of this script's folder; the scorer's block is read from <loop>/output/run.log
# (references/harness.md in the autoresearch skill fixes the block's keys).
#
# usage (any directory inside the repository):
#   bash <loop>/harness/record_result.sh <keep|discard|gate_failed|crash|timeout> "<description>" [--tie] [--killed]
#
# What it does, in order:
#   1. reads the block at the end of run.log and checks it belongs to HEAD (its commit and trial_md5 lines, which every
#      block carries, failed ones included) and that the trial file has no uncommitted change;
#   2. checks the status against the block: keep/discard need `status: ok`; gate_failed needs `status: gate_failed`;
#      crash/timeout need a failed block, or no block with --killed for a scoring the agent had to kill;
#   3. enforces the keep rule in the score's direction: gain = delta_vs_best (higher) or -delta_vs_best (lower);
#      keep needs gain > 0 and gain >= paired_se, or with --tie: 0 <= gain < paired_se (the trial file got simpler or
#      faster); a commit whose truth_hash differs from best_block.txt's is a truth change and needs gain >= 2 paired_se,
#      never a tie; the baseline (first row) is exempt and is never reset;
#   4. appends the row to <loop>/output/results.tsv (created with its header on the first call);
#   5. on keep copies last_scores.csv to best_scores.csv and the block to best_block.txt; otherwise runs
#      `git reset --keep HEAD~1` (uncommitted edits to other tracked files survive; output/ is untracked and survives);
#   6. prints the row, the git tip, and the running counts the stopping rules use.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$(git rev-parse --show-toplevel)"
top="$(pwd -P)"
loop_abs="$(dirname "$script_dir")"
loop_dir="${loop_abs#"$top"/}"   # the loop folder, relative to the repository root
conf="$script_dir/loop.conf"
[[ -f "$conf" ]] || { echo "$conf not found; write it (TRIAL_FILE, DIRECTION, EXTRA_FIELDS, SCORER)" >&2; exit 1; }
TRIAL_FILE=""; DIRECTION=""; EXTRA_FIELDS=""; SCORER="the scorer"
# shellcheck disable=SC1090
source "$conf"
[[ -n "$TRIAL_FILE" && -f "$TRIAL_FILE" ]] || { echo "loop.conf's TRIAL_FILE ('$TRIAL_FILE') is not a file at the repository root" >&2; exit 1; }
git ls-files --error-unmatch "$TRIAL_FILE" >/dev/null 2>&1 || { echo "$TRIAL_FILE is not tracked by git; commit the loop folder before logging" >&2; exit 1; }
case "$DIRECTION" in higher|lower) ;; *) echo "loop.conf's DIRECTION must be higher or lower (got '$DIRECTION')" >&2; exit 1 ;; esac

status="${1:-}"; description="${2:-}"; shift 2 || true
tie=0; killed=0
for flag in "$@"; do
  case "$flag" in
    --tie) tie=1 ;;
    --killed) killed=1 ;;
    *) echo "unknown flag $flag" >&2; exit 1 ;;
  esac
done
out="$loop_dir/output"
log="$out/run.log"
tsv="$out/results.tsv"
extra_header=""
for k in $EXTRA_FIELDS; do extra_header+=$'\t'"$k"; done
header=$'commit\ttime\tscore\tscore_se'"$extra_header"$'\tseconds\tstatus\tdescription'

case "$status" in
  keep|discard|gate_failed|crash|timeout) ;;
  *) echo "status must be keep, discard, gate_failed, crash or timeout (got '$status')" >&2; exit 1 ;;
esac
[[ -n "$description" ]] || { echo "a description is required" >&2; exit 1; }
[[ -f "$log" ]] || { echo "$log not found; run $SCORER dev first" >&2; exit 1; }
if [[ -f "$tsv" ]]; then
  [[ "$(head -n 1 "$tsv")" == "$header" ]] || { echo "$tsv has an unexpected header; fix it or loop.conf's EXTRA_FIELDS before logging" >&2; exit 1; }
fi
num_rows=0; [[ -f "$tsv" ]] && num_rows=$(($(wc -l < "$tsv") - 1))
is_baseline=0; [[ $num_rows -eq 0 ]] && is_baseline=1

# the block: everything after the last '---' line
block=$(awk '/^---$/{b=""; next} {b=b $0 "\n"} END{printf "%s", b}' "$log")
field() { printf '%s\n' "$block" | grep -m1 "^$1:" | awk '{print $2}' || true; }
head_commit=$(git rev-parse --short=7 HEAD)
log_status=$(field status); log_commit=$(field commit); log_md5=$(field trial_md5)
file_md5=$(md5 -q "$TRIAL_FILE" 2>/dev/null || md5sum "$TRIAL_FILE" | awk '{print $1}')

# the log must come from the committed trial file of HEAD; every block, failed ones included, carries both lines
if [[ -n "$log_status" && ( -z "$log_commit" || -z "$log_md5" ) ]]; then
  echo "run.log's block has no commit or trial_md5 line; the scorer prints both on every block (references/harness.md), rerun $SCORER" >&2; exit 1
fi
if [[ -n "$log_commit" && "$log_commit" != "$head_commit" ]]; then
  echo "run.log was produced at commit $log_commit but HEAD is $head_commit; rerun $SCORER on HEAD" >&2; exit 1
fi
if [[ -n "$log_md5" && "$log_md5" != "$file_md5" ]]; then
  echo "run.log was produced from a $TRIAL_FILE whose checksum differs from the file on disk; commit or revert, then rerun" >&2; exit 1
fi
if ! git diff --quiet HEAD -- "$TRIAL_FILE"; then
  echo "$TRIAL_FILE has uncommitted changes; amend them into the experiment commit or revert them first" >&2; exit 1
fi
if [[ $is_baseline -eq 0 ]] && cut -f1 "$tsv" | tail -n +2 | grep -qx "$head_commit"; then
  echo "commit $head_commit is already logged; a rescored experiment needs an amended commit (new hash)" >&2; exit 1
fi

# the status must match the block
case "$status" in
  keep|discard)
    [[ "$log_status" == "ok" ]] || { echo "run.log's status is '${log_status:-missing}', so log that status, not $status" >&2; exit 1; } ;;
  gate_failed)
    [[ "$log_status" == "gate_failed" ]] || { echo "run.log's status is '${log_status:-missing}', not gate_failed" >&2; exit 1; }
    [[ $is_baseline -eq 0 ]] || { echo "the baseline fails the gate: fix the baseline definition and rerun; nothing logged" >&2; exit 1; } ;;
  crash|timeout)
    if [[ "$log_status" == "ok" || "$log_status" == "gate_failed" ]]; then echo "run.log says status: $log_status; log keep, discard or gate_failed" >&2; exit 1; fi
    if [[ -z "$log_status" && $killed -eq 0 ]]; then
      echo "run.log has no block; pass --killed if you killed the scoring, otherwise it is still running" >&2; exit 1
    fi
    [[ $is_baseline -eq 0 ]] || { echo "the baseline failed: fix the bank or the environment and rerun; nothing logged" >&2; exit 1; } ;;
esac

# the keep rule
if [[ "$status" == "keep" && $is_baseline -eq 0 ]]; then
  delta=$(field delta_vs_best); pse=$(field paired_se)
  [[ "$delta" != "NA" && -n "$delta" ]] || { echo "delta_vs_best is NA: best_scores.csv is missing or from another bank; cannot keep" >&2; exit 1; }
  truth_now=$(field truth_hash); truth_best=""
  [[ -f "$out/best_block.txt" ]] && truth_best=$(grep -m1 '^truth_hash:' "$out/best_block.txt" | awk '{print $2}' || true)
  factor=1
  if [[ -n "$truth_now" && -n "$truth_best" && "$truth_now" != "$truth_best" ]]; then
    factor=2
    [[ $tie -eq 0 ]] || { echo "a truth change (truth_hash differs from the best commit's) is never kept as a tie" >&2; exit 1; }
  fi
  verdict=$(awk -v d="$delta" -v s="$pse" -v tie="$tie" -v dir="$DIRECTION" -v k="$factor" 'BEGIN{
    g = (dir == "higher") ? d + 0 : -(d + 0); s = s + 0;
    if (g > 0 && g >= k*s) print "keep";
    else if (tie == 1 && g >= 0 && (g < s || g == 0)) print "keep";
    else print "discard" }')
  if [[ "$verdict" != "keep" ]]; then
    echo "keep rule fails: delta_vs_best $delta, paired_se $pse, direction $DIRECTION, factor $factor (keep needs a gain of at least $factor paired_se, or --tie with a gain in [0, paired_se))" >&2
    exit 1
  fi
fi

# the row
extra_values=""
if [[ "$status" == "keep" || "$status" == "discard" || "$status" == "gate_failed" ]]; then
  score=$(field score); score_se=$(field score_se); seconds=$(field seconds)
  for v in "$score" "$score_se" "$seconds"; do
    [[ -n "$v" ]] || { echo "run.log's block is incomplete (score, score_se, seconds); rerun $SCORER" >&2; exit 1; }
  done
  for k in $EXTRA_FIELDS; do
    v=$(field "$k"); [[ -n "$v" ]] || { echo "run.log's block has no '$k:' line, which loop.conf's EXTRA_FIELDS names" >&2; exit 1; }
    extra_values+=$'\t'"$v"
  done
else
  score=0.000000; score_se=0.000000
  seconds=$(field seconds); seconds=${seconds:-0.0}
  for k in $EXTRA_FIELDS; do extra_values+=$'\t0'; done
fi
description=$(printf '%s' "$description" | tr '\t\n\r' '   ')
[[ -f "$tsv" ]] || printf '%s\n' "$header" > "$tsv"
printf '%s\t%s\t%s\t%s%s\t%s\t%s\t%s\n' \
  "$head_commit" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$score" "$score_se" "$extra_values" "$seconds" "$status" "$description" >> "$tsv"

# move the branch
if [[ "$status" == "keep" ]]; then
  cp "$out/last_scores.csv" "$out/best_scores.csv"
  printf '%s' "$block" > "$out/best_block.txt"
  echo "kept: best_scores.csv and best_block.txt now hold $head_commit"
elif [[ $is_baseline -eq 0 ]]; then
  git reset --keep HEAD~1 >/dev/null
  echo "reset: $TRIAL_FILE is back at $(git rev-parse --short=7 HEAD)"
fi

# the running counts the stopping rules use (experiments exclude the baseline row)
tail -n 1 "$tsv"
experiments=$(($(wc -l < "$tsv") - 2))
status_col=$(head -n 1 "$tsv" | tr '\t' '\n' | grep -nx status | cut -d: -f1)
streak=$(tail -n +3 "$tsv" | cut -f"$status_col" | awk '{ if ($1 == "keep") n = 0; else n++ } END { print n + 0 }')
first_time=$(sed -n '2p' "$tsv" | cut -f2)
echo "experiments after the baseline: $experiments; consecutive non-keeps: $streak; baseline logged at $first_time; tip $(git rev-parse --short=7 HEAD)"
