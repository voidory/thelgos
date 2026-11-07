import argparse,re
from pathlib import Path
import numpy as np
from bs4 import BeautifulSoup
from molsym.symtext.general_irrep_mats import pg_to_symels
from molsym.symtext.multiplication_table import build_mult_table
from common import ROOT,write_json,sha
JOBS={'Oh':(221,32),'Td':(215,31),'D3h':(187,26),'C2v':(25,7),'D2h':(47,8),'D4h':(123,15),'D2d':(111,14),'C2h':(10,5),'D6h':(191,27),'C3v':(156,19),'D8h':None}
I=np.eye(3);X=np.array([1,0,0]);Y=np.array([0,1,0]);Z=np.array([0,0,1])
def rot(u,t):
 u=np.array(u,float);u/=np.linalg.norm(u);t=np.deg2rad(t);x,y,z=u;A=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
 return np.cos(t)*I+(1-np.cos(t))*np.outer(u,u)+np.sin(t)*A
def mirror(u):
 u=np.array(u,float);u/=np.linalg.norm(u);return I-2*np.outer(u,u)
def norm(s,g,repair=True):
 s=s.replace('_','');s=re.sub(r"([AB])(''{0,1})([12])",r'\1\3\2',s)
 if g=='D2h':
  s=s.replace('A1','A')
  if repair:s={'B2g':'B3g','B3g':'B2g','B2u':'B3u','B3u':'B2u'}.get(s,s)
 return s
TARGETS={
'Oh':{'1':I,'4':rot(Z,90),'2100':rot(X,180),'3':rot([1,1,1],120),'2110':rot([1,1,0],180),'-1':-I,'-4':-rot(Z,90),'m100':mirror(X),'-3':-rot([1,1,1],120),'m110':mirror([1,1,0])},
'Td':{'1':I,'3':rot([1,1,1],120),'2100':rot(X,180),'-4':-rot(Z,90),'m110':mirror([1,1,0])},
'D3h':{'1':I,'m':mirror(Z),'3':rot(Z,120),'-6':-rot(Z,60),'2120':rot(X,180),'m100':mirror(Y)},
'C2v':{'1':I,'2z':rot(Z,180),'my':mirror(Y),'mx':mirror(X)},
# [K21] C3v (CONVENTIONS v1.5). Bilbao's m1-10 is in hexagonal axes: direction [1,-1,0] = Cartesian (3/2,-sqrt3/2,0); its mirror lies in the class of sigma_v(0), the xz plane.
'C3v':{'1':I,'3':rot(Z,120),'m1-10':mirror([np.sqrt(3)/2,-0.5,0])},
'D2h':{'1':I,'2z':rot(Z,180),'2y':rot(Y,180),'2x':rot(X,180),'-1':-I,'mz':mirror(Z),'my':mirror(Y),'mx':mirror(X)},
'D4h':{'1':I,'2':rot(Z,180),'4':rot(Z,90),'2100':rot(X,180),'21-10':rot([1,-1,0],180),'-1':-I,'mz':mirror(Z),'-4':-rot(Z,90),'m100':mirror(X),'m1-10':mirror([1,-1,0])},
'D2d':{'1':I,'2':rot(Z,180),'-4':-rot(Z,90),'2100':rot(X,180),'m1-10':mirror([1,-1,0])},
'C2h':{'1':I,'2':rot(Z,180),'-1':-I,'m':mirror(Z)},
'D6h':{'1':I,'6':rot(Z,60),'3':rot(Z,120),'2z':rot(Z,180),'2120':rot(X,180),'2100':rot(Y,180),'-1':-I,'-6':-rot(Z,60),'-3':-rot(Z,120),'mz':mirror(Z),'m120':mirror(X),'m100':mirror(Y)}}
EXACT=np.array(sorted(set([0.,1.,-1.,.5,-.5,2.,-2.,3.,-3.,np.sqrt(2),-np.sqrt(2),np.sqrt(2)/2,-np.sqrt(2)/2,np.sqrt(3)/2,-np.sqrt(3)/2])))
def snap(a):
 a=np.asarray(a,float);b=EXACT[np.argmin(abs(a[...,None]-EXACT),axis=-1)];e=float(np.max(abs(a-b)));assert e<1e-12,e;return b,e

def parse_bilbao(path,g):
 soup=BeautifulSoup(path.read_bytes(),'html.parser')
 for tab in soup.find_all('table'):
  rows=tab.find_all('tr',recursive=False) or tab.find_all('tr')
  if any('functions' in r.get_text().lower() for r in rows[:2]):break
 else:raise ValueError(g)
 rows=[[c.get_text('',strip=True).replace('−','-') for c in r.find_all(['td','th'],recursive=False)] for r in rows];heads=rows[0][2:-1];sizes=[1]*len(heads);vals={}
 for row in rows[1:]:
  if row[0]=='Mult.':sizes=[int(x) for x in row[2:-1]]
  elif len(row)==len(heads)+3:vals[norm(row[0],g,False)]=[float(x) for x in row[2:-1]]
 return heads,sizes,vals

def axis_label(mats,chars,k):
 v=mats[:,k,k];matches=[label for label,c in chars.items() if np.max(abs(c-v))<1e-12];assert len(matches)==1,matches;return matches[0]
def run(out,evidence=None):
 evidence=evidence or ROOT/'evidence'
 reports=[]
 for g,spec in JOBS.items():
  sym,irr,reps=pg_to_symels(g);M,me=snap([s.rrep for s in sym]);N=len(sym);T=build_mult_table(sym).astype(int)
  raw={norm(i.symbol,g,False):snap(np.trace(reps[i.symbol],axis1=1,axis2=2))[0] for i in irr}
  chars={norm(i.symbol,g):snap(np.trace(reps[i.symbol],axis1=1,axis2=2))[0] for i in irr}
  dims={norm(i.symbol,g):int(i.d) for i in irr};Xc=np.array(list(chars.values()));orth=float(np.max(abs(Xc@Xc.T/N-np.eye(len(irr)))))
  closure=max(float(np.max(abs(M[i]@M[j]-M[T[i,j]]))) for i in range(N) for j in range(N));assert closure<1e-12 and orth<1e-12
  ident=int(np.argmin(np.max(abs(M-I),axis=(1,2))));inv=[int(np.flatnonzero(T[i]==ident)[0]) for i in range(N)];remaining=set(range(N));classes=[]
  while remaining:
   i=min(remaining);cl=sorted(set(int(T[T[j,i],inv[j]]) for j in range(N)));classes.append(cl);remaining-=set(cl)
  gens=[];known={ident}
  for i in range(N):
   if i in known:continue
   gens.append(i)
   while True:
    new=known|{int(T[a,b]) for a in known for b in gens}
    if new==known:break
    known=new
  report={'group':g,'order':N,'closure_max_error':closure,'orthogonality_max_error':orth,'snap_matrix_max_residual':me,'snap_character_max_residual':max(snap(np.trace(reps[i.symbol],axis1=1,axis2=2))[1] for i in irr)}
  data={'schema_version':'1.2','group':g,'generator':'MolSym 1.3.1','operations':[{'symbol':s.symbol,'matrix':M[k].tolist()} for k,s in enumerate(sym)],'generator_indices':gens,'multiplication_table':T.tolist(),'classes':classes,'irreps':[{'label':label,'dimension':dims[label],'characters':c.tolist()} for label,c in chars.items()],'library_label_map':{norm(i.symbol,g,False):norm(i.symbol,g) for i in irr}}
  if spec:
   path=evidence/f'bilbao_{g}.html';heads,sizes,vals=parse_bilbao(path,g);mapping=[];idx=[]
   for h,n in zip(heads,sizes):
    target=TARGETS[g][h];hit=[cl for cl in classes if any(np.max(abs(M[j]-target))<1e-12 for j in cl)];assert len(hit)==1,(g,h,hit)
    cl=hit[0];assert len(cl)==n;idx.append(cl[0]);mapping.append({'bilbao_column':h,'target_matrix':snap(target)[0].tolist(),'operation_indices':cl,'representative_symbol':sym[cl[0]].symbol})
   assert len({tuple(x['operation_indices']) for x in mapping})==len(classes)
   assert set(vals)==set(chars)
   before=max(float(np.max(abs(raw[l][idx]-vals[l]))) for l in vals)
   err=max(float(np.max(abs(chars[l][idx]-vals[l]))) for l in vals);assert err<1e-12,(g,err)
   report.update(bilbao_status='PASS_GEOMETRY_MAPPED',bilbao_max_error=err,raw_library_max_error=before)
   sg,num=spec;data['bilbao']={'source':f'https://cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg={sg}&num={num}','snapshot_sha256':sha(path),'class_mapping':mapping}
   write_json(out/'checks'/f'bilbao_{g}.json',{'group':g,'class_mapping':mapping,'source_characters':vals,'comparison':report})
  else:report['bilbao_status']='NOT_COVERED_NONCRYSTALLOGRAPHIC'
  if g in ['D2h','C2v']:
   axis={a:axis_label(M,chars,k) for k,a in enumerate('xyz')};assert axis==({'x':'B3u','y':'B2u','z':'B1u'} if g=='D2h' else {'x':'B1','y':'B2','z':'A1'});report['axis_assertions']=axis
  if g in ['D3h','D6h']:
   xrot=rot([1,0,0],180);assert any(np.max(abs(m-xrot))<1e-12 for m in M)
   report['axis_assertions']={'C2prime_contains_x_axis':True}
   if g=='D6h':
    atom=np.array([1,0,0]);mid=np.array([np.cos(np.pi/6),np.sin(np.pi/6),0]);assert np.max(abs(xrot@atom-atom))<1e-12
    assert any(np.max(abs(m-rot(mid,180)))<1e-12 for m in M);report['benzene_atoms_on_pm_x']=True
  write_json(out/'groups'/f'{g}.json',data);reports.append(report)
 write_json(out/'checks'/'group_validation.json',reports)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--evidence',type=Path);a=p.parse_args();run(a.out,a.evidence)
