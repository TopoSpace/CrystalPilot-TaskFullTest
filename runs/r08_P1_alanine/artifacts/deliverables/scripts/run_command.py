"""Run one crystallographic command in a local subdirectory, recording provenance.
Usage: python run_command.py WORK_SUBDIR LABEL executable [arguments...]
No command is run in the read-only inputs directory.
"""
import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent

def run(work, label, command):
    work = (ROOT / work).resolve()
    if not work.is_relative_to(ROOT) or work.is_relative_to(ROOT / 'inputs'):
        raise ValueError('Work directory must be inside this run and outside inputs')
    work.mkdir(parents=True, exist_ok=True)
    log = work / (label + '.console.log')
    if log.exists():
        raise FileExistsError(log)
    env = os.environ.copy()
    prefix = Path(r'C:\users\<user>\miniforge3\envs\dials')
    env['PATH'] = str(prefix / 'Library' / 'bin') + os.pathsep + str(prefix / 'Scripts') + os.pathsep + env['PATH']
    env['PYTHONIOENCODING'] = 'utf-8'
    env['OMP_NUM_THREADS'] = '4'
    command[0] = shutil.which(command[0], path=env['PATH']) or command[0]
    record = {'start': datetime.datetime.now().isoformat(), 'cwd': str(work), 'argv': command, 'console_log': str(log)}
    start = time.monotonic()
    with log.open('w', encoding='utf-8') as f:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')
        f.flush()
        result = subprocess.run(command, cwd=work, env=env, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    record.update(exit_code=result.returncode, elapsed_seconds=round(time.monotonic()-start, 3), end=datetime.datetime.now().isoformat())
    with (ROOT / 'commands.jsonl').open('a', encoding='utf-8') as f:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')
    print(json.dumps(record, ensure_ascii=False, indent=2))
    lines = log.read_text(encoding='utf-8', errors='replace').splitlines()
    print('\n'.join(lines[-100:]))
    return result.returncode

if __name__ == '__main__':
    sys.exit(run(sys.argv[1], sys.argv[2], sys.argv[3:]))
