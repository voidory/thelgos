import os,sys,tempfile,subprocess,hashlib
from pathlib import Path
sys.dont_write_bytecode=True
from common import ROOT

def tree_hash(root):return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
def main():
 before=tree_hash(ROOT);env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OPENBLAS_CORETYPE':'Haswell'}  # [K20] same kernel on every host
 with tempfile.TemporaryDirectory(prefix='lgo12-baseline-') as tmp:
  out=Path(tmp)
  for script in ['generate_groups.py','check_release1.py','build_checks.py']:
   subprocess.run([sys.executable,str(ROOT/'scripts'/script),'--out',str(out)],env=env,check=True)
  # Compare exact serialized outputs; no committed data is overwritten.
  for name in ['groups','checks']:
   subprocess.run(['diff','-r',str(ROOT/name),str(out/name)],check=True)
  subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=ROOT,check=True)
 assert before==tree_hash(ROOT),'check modified committed bundle'
 print('PASS: temporary regeneration, diff -r groups/checks, sha256sum and no mutation')
if __name__=='__main__':main()
