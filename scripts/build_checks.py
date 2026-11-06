import argparse
from pathlib import Path
from common import ROOT,write_json,sha
from check_support import run
from compare_crystod import run as crystod

def build(out):
 run(out)
 code=crystod(out,cases_file=out/'checks/crystod_cases.json');assert code==2,'unavailable CrystOD must block, not pass'
 import json
 lock=json.loads((ROOT/'evidence/pyodide_lock.txt').read_text());assert lock['info']['python']=='3.12.7' and lock['info']['version']=='0.27.7' and lock['packages']['numpy']['version']=='2.0.2'
 write_json(out/'checks/pyodide_pin.json',{'status':'PASS_METADATA_ONLY','checked_on':'2026-11-07','python':lock['info']['python'],'numpy':lock['packages']['numpy']['version'],'pyodide':lock['info']['version'],'source':'https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide-lock.json','package_page':'https://pyodide.org/en/0.27.7/usage/packages-in-pyodide.html','snapshot_sha256':sha(ROOT/'evidence/pyodide_lock.txt'),'browser_execution':'NOT_RUN'})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);build(p.parse_args().out)
