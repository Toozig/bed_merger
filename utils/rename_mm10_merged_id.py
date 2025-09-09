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



def main():
    bed_df = pd.read_csv(BED_PATH, sep='\t', header=None)
    xl_df = pd.read_excel(EXCEL)
    xl_df['new_id'] = 'mm10_' + xl_df['id'].apply(shorten_developmental_id) + '.' + xl_df.index.astype(str)
    bed_df[3] = xl_df['new_id']
    bed_df.to_csv(BED_PATH, sep='\t', header=False, index=False)

if __name__ == '__main__':
    main()