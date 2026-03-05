#!/usr/bin/env bash

set -euo pipefail

# This script scans the project's results/ directory for dataset folders,
# excluding Supp_Data_3*, old_result, and human_yaml_run. For each dataset,
# it finds the first *Statistics_Summary.csv, extracts the first data row's
# Display Name (second line), generates pie charts, and uploads the figures
# to Dropbox under the requested destination path.

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
RESULTS_DIR="$ROOT_DIR/results"

DEST_ROOT="Nitzan_Gonen_lab/Joint_projects/Lab_Manuscripts/WGS_DSD_paper/Tables"

# Embed upload_dir_to_dropbox to avoid requiring external sourcing
upload_dir_to_dropbox() {
    local source_dir="$1"
    local dest_path="$2"

    if [[ $# -ne 2 ]]; then
        echo "Usage: upload_dir_to_dropbox <source_directory> <destination_path>"
        echo "Example: upload_dir_to_dropbox results/figures/ 'Nitzan_Gonen_lab/project/figures/'"
        return 1
    fi
    if [[ ! -d "$source_dir" ]]; then
        echo "Error: Source directory '$source_dir' does not exist"
        return 1
    fi
    if ! command -v dbxcli &> /dev/null; then
        echo "Error: dbxcli command not found. Please install Dropbox CLI."
        return 1
    fi
    echo "Uploading contents of '$source_dir' to '$dest_path'..."
    for file in "$source_dir"/*; do
        if [[ -f "$file" ]]; then
            local filename
            filename=$(basename "$file")
            echo "Uploading: $filename"
            dbxcli put "$file" "$dest_path/$filename"
        fi
    done
    echo "Upload completed."
}

# Use project virtual environment if present
if [[ -d "$ROOT_DIR/.venv" ]]; then
  # shellcheck disable=SC1091
  source "$ROOT_DIR/.venv/bin/activate"
fi

shopt -s nullglob

for dataset_dir in "$RESULTS_DIR"/*/; do
  base="$(basename "$dataset_dir")"

  # Exclusions
  if [[ "$base" == old_result ]] || [[ "$base" == human_yaml_run ]] || [[ "$base" == Supp_Data_3* ]]; then
    continue
  fi

  # Find the first Statistics Summary CSV in this dataset directory
  csv_candidates=( "$dataset_dir"/*Statistics_Summary.csv )
  if (( ${#csv_candidates[@]} == 0 )); then
    continue
  fi
  csv="${csv_candidates[0]}"

  # Extract the Display Name from the first data row (second line)
  # Trim surrounding quotes and whitespace
  display_name=$(awk -F',' 'NR==2 {print $1; exit}' "$csv" | sed -e 's/^"//' -e 's/"$//' -e 's/^\s*//' -e 's/\s*$//')

  if [[ -z "$display_name" ]]; then
    # If empty, skip this dataset
    continue
  fi

  # Generate figures into <dataset_dir>/figures using the Python CLI
  python "$ROOT_DIR/utils/plot_accessibility_pie.py" "$csv" "$dataset_dir" --rows "$display_name"

  # Upload the figures directory to Dropbox under <DEST_ROOT>/<dataset>/figures
  src="$dataset_dir/figures"
  if [[ -d "$src" ]]; then
    dest="$DEST_ROOT/$base/figures"
    upload_dir_to_dropbox "$src" "$dest"
  fi
done

echo "Done."


