import csv,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def write_json(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def csv_read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
