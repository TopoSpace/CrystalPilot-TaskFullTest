from pathlib import Path
import json
from merge import merge

OUT=Path(__file__).resolve().parent
merge("unmerged_absorption.csv","scale_validation_constant",validate_scales=True)
merge("unmerged_absorption.csv","scale_validation_scan",scan=True,validate_scales=True)
a=json.loads((OUT/"scale_validation_constant_statistics.json").read_text())
b=json.loads((OUT/"scale_validation_scan_statistics.json").read_text())
scan=b["scale_validation_Rint_filtered"] < .97*a["scale_validation_Rint_filtered"]
print("chosen scan polynomial",scan,flush=True)
m,df=merge("unmerged_absorption.csv","merged_corrected",scan=scan)
selected=m[(~m.absent_212121)&(m.d>=.75)]
selected.to_csv(OUT/"merged_final.csv",index=False)
selected.to_csv(OUT/"phasing_input.csv",index=False)
