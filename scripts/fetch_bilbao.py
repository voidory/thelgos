import argparse
from pathlib import Path
from urllib.request import urlopen
from generate_groups import JOBS
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 for group,spec in JOBS.items():
  if not spec:continue
  sg,num=spec;url=f'https://cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg={sg}&num={num}'
  with urlopen(url,timeout=60) as r:content=r.read()
  assert b'Character Table' in content;(a.out/f'bilbao_{group}.html').write_bytes(content)
 print('Fresh Bilbao snapshots fetched; no committed file overwritten unless explicitly used as --out')
