#!/usr/bin/env bash
set -euo pipefail

# Run coercive-field calculations for all FE2 switching paths in one phase folder.
#
# Expected layout:
#   fe2_phase_folder/
#     xxx-nnn-nnn/
#       00/ 01/ ... 1-dipol/ a-xxx-nnn-nnn/ b-xxx-nnn-nnn/ c-xxx-nnn-nnn/
#     yxz-nnn-nnn/
#     compute_coercive_field/        # this script lives here
#
# Path names are matched by: [xyz][xyz][xyz]-[01][01][01]-[01][01][01]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$SCRIPT_DIR"
PHASE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_ROOT="$PHASE_ROOT/coercive_field_all_results"
PYTHON_BIN="${PYTHON_BIN:-python}"
EXTRA_ARGS_TEXT=""
DRY_RUN=0
OUTPUT_ROOT_EXPLICIT=0

usage() {
    cat <<EOF
Usage:
  bash run_all_fe2_paths.sh [options]

Options:
  --phase-root DIR       Folder containing all xxx-nnn-nnn paths.
                         Default: parent folder of this script.
  --code-dir DIR         Folder containing compute_coercive_field.py.
                         Default: this script folder.
  --output-root DIR      Folder for all outputs.
                         Default: <phase-root>/coercive_field_all_results.
  --python CMD           Python command. Default: python or \$PYTHON_BIN.
  --extra-args "ARGS"    Extra arguments passed to compute_coercive_field.py.
                         Example: --extra-args "--initial-polarization 0,0,50"
  --dry-run              Print commands without running them.
  -h, --help             Show this help.

Examples:
  bash run_all_fe2_paths.sh
  bash run_all_fe2_paths.sh --python python3
  bash run_all_fe2_paths.sh --phase-root /path/to/fe2 --code-dir /path/to/fe2/comput_coercive_field
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --phase-root)
            PHASE_ROOT="$(cd "$2" && pwd)"
            shift 2
            ;;
        --code-dir)
            CODE_DIR="$(cd "$2" && pwd)"
            shift 2
            ;;
        --output-root)
            mkdir -p "$2"
            OUTPUT_ROOT="$(cd "$2" && pwd)"
            OUTPUT_ROOT_EXPLICIT=1
            shift 2
            ;;
        --python)
            PYTHON_BIN="$2"
            shift 2
            ;;
        --extra-args)
            # shellcheck disable=SC2206
            EXTRA_ARGS_TEXT="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ "$OUTPUT_ROOT_EXPLICIT" -eq 0 ]]; then
    OUTPUT_ROOT="$PHASE_ROOT/coercive_field_all_results"
fi

if [[ ! -f "$CODE_DIR/compute_coercive_field.py" ]]; then
    echo "Cannot find compute_coercive_field.py in code dir: $CODE_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_ROOT"
LOG_FILE="$OUTPUT_ROOT/run_all_fe2_paths.log"
SUMMARY_FILE="$OUTPUT_ROOT/all_coercive_field_summary.csv"
STATUS_FILE="$OUTPUT_ROOT/run_status.tsv"

: > "$LOG_FILE"
printf "path\tfield_tag\tfield_direction\tstatus\toutput_dir\n" > "$STATUS_FILE"
rm -f "$SUMMARY_FILE"

# Common field directions:
#   001:  Pca2_1 ferroelectric axis, using the field opposite to initial [001].
#   010:  common non-180 Pca2_1 endpoint direction.
#   101:  Pmn2_1 ferroelectric axis, using the field opposite to initial [101].
#   m110: common non-180 Pmn2_1 endpoint [-110], i.e. (-1,1,0).
#   111:  retained for comparison with older direction scans.
FIELD_DIRECTIONS=(
    "001:0,0,-1"
    "010:0,1,0"
    "101:-1,0,-1"
    "m110:-1,1,0"
    "111:-1,-1,-1"
)

echo "Code dir    : $CODE_DIR" | tee -a "$LOG_FILE"
echo "Phase root  : $PHASE_ROOT" | tee -a "$LOG_FILE"
echo "Output root : $OUTPUT_ROOT" | tee -a "$LOG_FILE"
echo "Python      : $PYTHON_BIN" | tee -a "$LOG_FILE"

mapfile -t PATH_DIRS < <(
    find "$PHASE_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\t%p\n' \
        | awk -F '\t' '$1 ~ /^[xyz][xyz][xyz]-[01][01][01]-[01][01][01]$/ {print $2}' \
        | sort
)

if [[ "${#PATH_DIRS[@]}" -eq 0 ]]; then
    echo "No FE2 path folders found under: $PHASE_ROOT" >&2
    echo "Expected names like xyz-101-000." >&2
    exit 1
fi

append_summary() {
    local path_name="$1"
    local field_tag="$2"
    local field_direction="$3"
    local summary_csv="$4"

    if [[ ! -f "$summary_csv" ]]; then
        return
    fi

    if [[ ! -f "$SUMMARY_FILE" ]]; then
        head -n 1 "$summary_csv" | awk -v OFS=',' '{print "path_name,field_tag,field_direction," $0}' > "$SUMMARY_FILE"
    fi

    tail -n +2 "$summary_csv" \
        | awk -v OFS=',' -v p="$path_name" -v t="$field_tag" -v d="$field_direction" '{print p,t,d,$0}' \
        >> "$SUMMARY_FILE"
}

run_one() {
    local path_dir="$1"
    local field_tag="$2"
    local field_direction="$3"
    local path_name
    local output_dir
    local summary_csv
    local extra_args=()

    path_name="$(basename "$path_dir")"
    output_dir="$OUTPUT_ROOT/${path_name}/field_${field_tag}"
    summary_csv="$output_dir/coercive_field_summary.csv"

    if [[ ! -d "$path_dir/00" || ! -d "$path_dir/1-dipol" ]]; then
        echo "Skip $path_name: missing 00/ or 1-dipol/" | tee -a "$LOG_FILE"
        printf "%s\t%s\t%s\tSKIPPED_MISSING_DATA\t%s\n" "$path_name" "$field_tag" "$field_direction" "$output_dir" >> "$STATUS_FILE"
        return
    fi

    mkdir -p "$output_dir"

    local cmd=(
        "$PYTHON_BIN" "$CODE_DIR/compute_coercive_field.py"
        --root "$path_dir"
        --field-direction="$field_direction"
        --output-dir "$output_dir"
    )

    if [[ -n "$EXTRA_ARGS_TEXT" ]]; then
        # shellcheck disable=SC2206
        extra_args=($EXTRA_ARGS_TEXT)
        cmd+=("${extra_args[@]}")
    fi

    echo "" | tee -a "$LOG_FILE"
    echo "Running path=$path_name field=$field_direction output=$output_dir" | tee -a "$LOG_FILE"
    printf 'Command:' | tee -a "$LOG_FILE"
    printf ' %q' "${cmd[@]}" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf "%s\t%s\t%s\tDRY_RUN\t%s\n" "$path_name" "$field_tag" "$field_direction" "$output_dir" >> "$STATUS_FILE"
        return
    fi

    if "${cmd[@]}" 2>&1 | tee -a "$LOG_FILE"; then
        append_summary "$path_name" "$field_tag" "$field_direction" "$summary_csv"
        printf "%s\t%s\t%s\tOK\t%s\n" "$path_name" "$field_tag" "$field_direction" "$output_dir" >> "$STATUS_FILE"
    else
        printf "%s\t%s\t%s\tFAILED\t%s\n" "$path_name" "$field_tag" "$field_direction" "$output_dir" >> "$STATUS_FILE"
        return 1
    fi
}

fail_count=0
for path_dir in "${PATH_DIRS[@]}"; do
    for item in "${FIELD_DIRECTIONS[@]}"; do
        field_tag="${item%%:*}"
        field_direction="${item#*:}"
        if ! run_one "$path_dir" "$field_tag" "$field_direction"; then
            fail_count=$((fail_count + 1))
        fi
    done
done

echo "" | tee -a "$LOG_FILE"
echo "Combined summary: $SUMMARY_FILE" | tee -a "$LOG_FILE"
echo "Run status      : $STATUS_FILE" | tee -a "$LOG_FILE"
echo "Log file        : $LOG_FILE" | tee -a "$LOG_FILE"

if [[ "$fail_count" -gt 0 ]]; then
    echo "Finished with $fail_count failed calculation(s)." >&2
    exit 1
fi

echo "Finished all calculations successfully." | tee -a "$LOG_FILE"
