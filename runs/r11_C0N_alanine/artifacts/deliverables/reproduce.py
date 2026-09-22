"""Run only the interpreter authorized in environment.json."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent


def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--full",action="store_true")
    mode.add_argument("--refine-only",action="store_true")
    mode.add_argument("--check-only",action="store_true")
    args=parser.parse_args()
    authorized=Path(json.loads((ROOT/"environment.json").read_text(encoding="utf-8"))["environment"]["python"])
    if Path(sys.executable).resolve()!=authorized.resolve():
        raise RuntimeError("Use only the Python interpreter in environment.json")
    if args.full:
        commands=[["decode_all.py"],["fit_geometry.py"],["geometry_flexible.py"],
                  ["integrate.py"],["shadow_mask.py"],["integrate_flexible.py"],
                  ["absorption.py"],["final_merge.py"],["solve.py"],["refine.py"],
                  ["refine_aniso.py","--final"],["validate.py"],["export_results.py"],
                  ["check_outputs.py"]]
    elif args.refine_only:
        commands=[["refine_aniso.py","--final"],["validate.py"],["export_results.py"],
                  ["check_outputs.py"]]
    else:
        commands=[["check_outputs.py"]]
    env=os.environ.copy()
    env["OPENBLAS_NUM_THREADS"]="1"
    env["OMP_NUM_THREADS"]="1"
    env["PYTHONHASHSEED"]="0"
    logs=OUT/"reproduction_logs"
    logs.mkdir(exist_ok=True)
    results=[]
    for command in commands:
        start=time.monotonic()
        print("RUN"," ".join(command),flush=True)
        with (logs/(Path(command[0]).stem+".log")).open("w",encoding="utf-8") as log:
            proc=subprocess.Popen([str(authorized),"-u",str(OUT/command[0]),*command[1:]],
                                  cwd=ROOT,env=env,stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT,text=True,encoding="utf-8",errors="replace")
            for line in proc.stdout:
                log.write(line)
                log.flush()
                print(line,end="",flush=True)
            code=proc.wait()
        results.append({"command":command,"returncode":code,"seconds":time.monotonic()-start})
        if code:
            raise RuntimeError(f"Stage failed: {command}")
    (OUT/"reproduction_run.json").write_text(json.dumps(results,indent=2))
    from check_outputs import provenance
    provenance()


if __name__=="__main__":
    main()
