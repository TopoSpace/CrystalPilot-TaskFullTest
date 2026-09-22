"""Reproduce this case with the explicitly authorized Python interpreter only."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parent
HASHES = {
    "start.hkl": "52e803152c644abbfc5a52d70ce90998b439efafec3922a15e4b7a461d8f2614",
    "start.ins": "f89c0cfbf199126ad800d146a6d1ef66cee594a537258245600056bad11207fb",
    "synthesis.png": "942b46860d0bf98be956c5ec33b731bcc03c836d37a470f95e8abaad31c0cc9f",
}

CORE_COMMANDS = [
    ["framework.py", "merge"],
    ["framework.py", "patterson"],
    ["solve_heavy.py"],
    ["inspect_maps.py"],
    ["refine.py"],
    ["refine.py", "--input", "isotropic_atoms.json", "--name", "aniso_heavy",
     "--aniso", "heavy", "--cycles", "5"],
    ["refine.py", "--input", "aniso_heavy_atoms.json", "--name", "aniso_all",
     "--aniso", "all", "--cycles", "5"],
    ["refine.py", "--input", "aniso_all_atoms.json", "--name", "highangle",
     "--aniso", "all", "--zr-fp", "--restraints", "--dmax", "3", "--cycles", "5"],
    ["refine.py", "--input", "highangle_atoms.json", "--name", "riding",
     "--aniso", "all", "--zr-fp", "--restraints", "--hydrogen", "--dmax", "3",
     "--cycles", "6"],
    ["audit_data.py"],
    ["refine.py", "--input", "riding_atoms.json", "--name", "final_all_data",
     "--dmin", "0.69", "--aniso", "all", "--zr-fp", "--restraints", "--hydrogen",
     "--cycles", "7"],
    ["refine.py", "--input", "riding_atoms.json", "--name", "highangle_check",
     "--dmin", "1.0", "--dmax", "3", "--aniso", "all", "--zr-fp", "--restraints",
     "--hydrogen", "--cycles", "6"],
    ["refine.py", "--input", "final_all_data_atoms.json", "--name", "no_dispersion",
     "--dmin", "0.69", "--aniso", "all", "--fix-zr-fp", "0", "--restraints",
     "--hydrogen", "--cycles", "6"],
    ["refine.py", "--input", "riding_atoms.json", "--name", "holdout_check",
     "--dmin", "0.69", "--aniso", "all", "--zr-fp", "--restraints", "--hydrogen",
     "--holdout", "--cycles", "6"],
]

EXPERIMENT_COMMANDS = [
    ["refine.py", "--input", "aniso_all_atoms.json", "--name", "aniso_fp",
     "--aniso", "all", "--zr-fp", "--cycles", "5"],
    ["refine.py", "--input", "aniso_all_atoms.json", "--name", "restrained",
     "--aniso", "all", "--restraints", "--cycles", "5"],
    ["solvent_density.py", "--input", "highangle", "--name", "pore_cv",
     "--holdout", "--iterations", "201"],
    ["bulk_pore.py", "--input", "highangle", "--kind", "layered", "--name", "bulk_layered"],
    ["bulk_pore.py", "--input", "riding", "--kind", "radial", "--name", "bulk_radial"],
    ["solvent_density.py", "--input", "restrained", "--name", "pore_conservative_cv",
     "--holdout", "--iterations", "81", "--cap", "0.5", "--smooth", "0.3",
     "--dmin", "3"],
]

END_COMMANDS = [
    ["diagnostics.py", "final_all_data"],
    ["diagnostics.py", "highangle_check"],
    ["export_results.py"],
    ["validate_exports.py"],
]


def verify_inputs():
    actual = {name:hashlib.sha256((WORKSPACE/"inputs"/name).read_bytes()).hexdigest()
              for name in HASHES}
    if actual != HASHES:
        raise RuntimeError("Input hashes differ from the supplied case")
    return actual


def compact_generated_arrays():
    import numpy as np
    model_names = ["isotropic", "aniso_heavy", "aniso_all", "aniso_fp", "restrained",
                   "highangle", "riding", "final_all_data", "highangle_check",
                   "no_dispersion", "holdout_check"]
    generated_maps = [f"{name}_{kind}.npy" for name in model_names
                      for kind in ("fo", "diff", "2fofc")]
    generated_maps += ["patterson.npy", "heavy_fourier.npy", "heavy_difference.npy"]
    compacted = []
    for filename in generated_maps:
        path = ROOT/filename
        if not path.exists():
            continue
        if path.resolve().parent != ROOT.resolve():
            raise RuntimeError("Output path left the deliverables directory")
        if filename in ("patterson.npy", "heavy_fourier.npy", "heavy_difference.npy"):
            density = np.load(path).astype(np.float32)
            np.savez_compressed(path.with_suffix(".npz"), density=density)
        # All final electron-density maps are retained together in final_maps.npz;
        # intermediate map maxima and their atom models remain in JSON files.
        path.unlink()
        compacted.append(filename)
    return compacted


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiments", action="store_true",
                        help="Also repeat rejected pore-density and scattering trials")
    parser.add_argument("--keep-large-maps", action="store_true")
    parser.add_argument("--compact-only", action="store_true")
    args = parser.parse_args()
    env = json.loads((WORKSPACE/"environment.json").read_text(encoding="utf-8-sig"))
    authorized = os.path.normcase(os.path.abspath(env["environment"]["python"]))
    if os.path.normcase(os.path.abspath(sys.executable)) != authorized:
        raise RuntimeError("Use only the Python path authorized by environment.json")
    before = verify_inputs()
    if args.compact_only:
        print(json.dumps(dict(compacted=compact_generated_arrays()), indent=2))
        return
    start = time.time()
    commands = CORE_COMMANDS + (EXPERIMENT_COMMANDS if args.experiments else []) + END_COMMANDS
    manifest = dict(start=datetime.datetime.now().astimezone().isoformat(),
                    interpreter=sys.executable, python_version=sys.version,
                    inputs=before, commands=[], external_downloads=False,
                    crystallographic_software_used=False)
    with (ROOT/"reproduction.log").open("w", encoding="utf-8") as log:
        for number, args_list in enumerate(commands, start=1):
            command = [sys.executable, str(ROOT/args_list[0]), *args_list[1:]]
            step_start = time.time()
            stamp = datetime.datetime.now().astimezone().isoformat()
            print(f"[{number}/{len(commands)}] {' '.join(args_list)}", flush=True)
            log.write(f"\n{stamp}\nCOMMAND {json.dumps(command)}\n")
            log.flush()
            process = subprocess.Popen(command, cwd=WORKSPACE, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True, encoding="utf-8")
            for line in process.stdout:
                log.write(line)
            code = process.wait()
            log.flush()
            manifest["commands"].append(dict(argv=command, exit_code=code,
                                               started=stamp, seconds=time.time()-step_start))
            (ROOT/"reproduction_manifest.json").write_text(json.dumps(manifest, indent=2))
            if code:
                raise RuntimeError(f"Calculation step failed: {args_list}")
    manifest["inputs_unchanged"] = verify_inputs() == before
    if not args.keep_large_maps:
        manifest["compacted_generated_maps"] = compact_generated_arrays()
    manifest["finished"] = datetime.datetime.now().astimezone().isoformat()
    manifest["seconds"] = time.time()-start
    (ROOT/"reproduction_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Complete: {len(commands)} stages, {time.time()-start:.1f} seconds. Inputs unchanged.")


if __name__ == "__main__":
    main()
