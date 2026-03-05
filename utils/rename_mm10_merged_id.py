import pandas as pd

from typing import List, Set
import re

import re
from typing import List, Set

BED_PATH = 'results/Supp_data_2_mm10_gnoad_cells_ATAC_seq/Supp_data_2_mm10_gnoad_cells_ATAC_seq.bed'
EXCEL = 'results/Supp_data_2_mm10_gnoad_cells_ATAC_seq/Supp_data_2_mm10_gnoad_cells_ATAC_seq.xlsx'

def shorten_developmental_id(long_id: str) -> str:
    """
    Convert long developmental IDs to a compact string.

    Supports both token orders:
    - E<stage>.<sex>.<number>
    - <sex>.E<stage>.<number>

    Example:
      "XX.E15.5.1_XX.E13.5.1" -> "XX13.5.15.5"

    Returns:
      "XY<sorted stages>XX<sorted stages>" (omit a group if absent)
    """
    pattern_sex_second = r"E(\d+\.?\d*)\.([XY]{2})\.\d+"
    pattern_sex_first = r"([XY]{2})\.E(\d+\.?\d*)\.\d+"

    # Collect (sex, stage) pairs from both formats
    matches: List[tuple[str, str]] = []
    matches.extend([(sex, stage) for stage, sex in re.findall(pattern_sex_second, long_id)])
    matches.extend([(sex, stage) for sex, stage in re.findall(pattern_sex_first, long_id)])

    if not matches:
        return long_id

    xy_stages: Set[str] = set()
    xx_stages: Set[str] = set()
    for sex, stage in matches:
        #remove 'E' & '.5'
        stage = stage.replace('E', '').replace('.5', '')
        if sex == "XY":
            xy_stages.add(stage)
        elif sex == "XX":
            xx_stages.add(stage)

    def sort_stages(stages: Set[str]) -> List[str]:
        return sorted(stages, key=lambda x: float(x))

    result_parts: List[str] = []
    if xy_stages:
        result_parts.append("XY" + ".".join(sort_stages(xy_stages)))
    if xx_stages:
        result_parts.append("XX" + ".".join(sort_stages(xx_stages)))

    return "".join(result_parts)


def insert_sheet_to_excel(
    excel_path: str,
    df: pd.DataFrame,
    sheet_name: str,
    overwrite: bool = True,
    index: bool = False
) -> None:
    """
    Insert a DataFrame into a specific sheet of an existing Excel file.
    
    Args:
        excel_path: Path to the Excel file
        df: DataFrame to insert
        sheet_name: Name of the sheet to create/overwrite
        overwrite: If True, overwrite existing sheet; if False, raise error if sheet exists
        index: Whether to include DataFrame index in the output
    
    Raises:
        FileNotFoundError: If Excel file doesn't exist
        ValueError: If sheet exists and overwrite=False
    """
    # Read existing Excel file to get all sheets
    existing_sheets = pd.ExcelFile(excel_path).sheet_names
    
    # Check if sheet already exists
    if sheet_name in existing_sheets and not overwrite:
        raise ValueError(f"Sheet '{sheet_name}' already exists and overwrite=False")
    
    # Create a copy of existing sheets data
    sheets_data = {}
    for sheet in existing_sheets:
        if sheet != sheet_name:  # Don't read the sheet we're about to overwrite
            sheets_data[sheet] = pd.read_excel(excel_path, sheet_name=sheet)
    
    # Add the new sheet data
    sheets_data[sheet_name] = df
    
    # Write all sheets back to Excel
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        for sheet, data in sheets_data.items():
            data.to_excel(writer, sheet_name=sheet, index=index)



def main():
    bed_df = pd.read_csv(BED_PATH, sep='\t', header=None)
    xl_df = pd.read_excel(EXCEL, sheet_name='merged')
    old_id = xl_df['id']
    new_id = 'mm10_' + xl_df['id'].apply(shorten_developmental_id) + '.' + xl_df.index.astype(str)
    xl_df['id'] = new_id
    xl_df['original_peaks'] = old_id
    bed_df[3] = new_id
    bed_df['original_peaks'] = old_id
    insert_sheet_to_excel(EXCEL, xl_df, 'merged', overwrite=True, index=False)
    xl_df.to_excel(EXCEL, sheet_name='merged', index=False)
    bed_df.to_csv(BED_PATH, sep='\t', header=False, index=False)

if __name__ == '__main__':
    main()