import argparse,csv,json
from pathlib import Path
import numpy as np
from common import ROOT,write_json,csv_read
from parameters import Parameters

def E2(a,b,g):return np.linalg.eigvalsh([[a,g],[g,b]]).tolist()
def warn_shell(d,n):
 partial=0<n<2*d
 return ['occupation_not_unique','Jahn_Teller_warning','spin_state_not_determined'] if d>1 and partial else []

def run(data,out,groups,check_reference=True):
 P=Parameters(data);reports=[];rows=[];refs={};cases=[]
 def test(name,group,atoms,basis,preset,ne,expected):
  p=P.presets[preset];N=len(basis);H=np.zeros((N,N));edges=[]
  for i,(a,v) in enumerate(basis):H[i,i]=P.onsite(atoms[a][0],'s' if v is None else 'p')
  for a in range(len(atoms)):
   for b in range(a+1,len(atoms)):
    u=atoms[b][1]-atoms[a][1];d=float(np.linalg.norm(u));u/=d
    if len(atoms)==2:typ='all'
    elif group=='C2v' and name=='allyl cation pi':typ='nearest_neighbor' if b-a==1 else 'absent'
    elif atoms[0][0] in ['O','B','C','S'] and name not in ['allyl cation pi']:typ='central_ligand' if a==0 else ('trans' if np.dot(atoms[a][1],atoms[b][1])<-.99*np.linalg.norm(atoms[a][1])*np.linalg.norm(atoms[b][1]) else 'cis')
    elif preset=='sixH_lecture':typ='trans' if np.dot(atoms[a][1],atoms[b][1])<-.99 else 'cis'
    else:typ='all'
    on=p['edge_policy']['edges_on'];off=p['edge_policy']['edges_off'];enabled=(typ in on or 'all' in on) and typ not in off
    if typ in ['cis','trans'] and 'ligand_ligand' in off:enabled=False
    if not enabled:continue
    edges.append([a,b]);block=P.block(atoms[a][0],atoms[b][0],u,d,preset)
    for i,(ia,v) in enumerate(basis):
     if ia!=a:continue
     va=np.array([1.]) if atoms[a][0]=='H' else (np.r_[1.,[0,0,0]] if v is None else np.r_[0.,v])
     for j,(ib,w) in enumerate(basis):
      if ib!=b:continue
      vb=np.array([1.]) if atoms[b][0]=='H' else (np.r_[1.,[0,0,0]] if w is None else np.r_[0.,w])
      H[i,j]=H[j,i]=va@block@vb
  e,c=np.linalg.eigh(H);analytic_err=float(np.max(abs(e-np.sort(expected))));assert analytic_err<1e-10,(name,analytic_err)
  g=json.loads((groups/f'{group}.json').read_text());ops=[]
  for op in g['operations']:
   m=np.array(op['matrix']);D=np.zeros((N,N))
   for j,(a,v) in enumerate(basis):
    dst=[b for b,(el,pos) in enumerate(atoms) if el==atoms[a][0] and np.linalg.norm(pos-m@atoms[a][1])<1e-9];assert len(dst)==1
    for i,(b,w) in enumerate(basis):
     if b==dst[0]:
      if v is None and w is None:D[i,j]=1
      elif v is not None and w is not None:D[i,j]=w@m@v
   assert np.max(abs(D.T@D-np.eye(N)))<1e-10;ops.append(D)
  ops=np.array(ops);projectors={r['label']:r['dimension']/len(ops)*np.einsum('g,gij->ij',r['characters'],ops) for r in g['irreps']}
  comm=float(max(np.max(abs(H@D-D@H)) for D in ops));assert comm<1e-10
  k=0;left=ne;levels=[];allwarnings=[]
  while k<N:
   ix=np.flatnonzero(abs(e-e[k])<1e-8);ix=ix[ix>=k];V=c[:,ix];weights={l:float(np.trace(V.T@pr@V)/len(ix)) for l,pr in projectors.items()};l=max(weights,key=weights.get);assert weights[l]>.999,(name,weights)
   d=len(ix);n=min(left,2*d);left-=n;warnings=warn_shell(d,n);allwarnings.extend(warnings)
   row=dict(case=name,irrep=l,energy_eV=float(e[k]),multiplicity=d,electrons=n,warnings=';'.join(warnings));rows.append(row);levels.append(row);k=int(ix[-1])+1
  assert left==0
  resid=float(np.linalg.norm(H@c-c*e,axis=0).max());ortho=float(np.max(abs(c.T@c-np.eye(N))));assert resid<1e-10 and ortho<1e-10
  refs[name]=[dict(irrep=x['irrep'],energy_eV=x['energy_eV'],multiplicity=x['multiplicity'],electrons=x['electrons']) for x in levels]
  report=dict(case=name,status='PASS',analytic_max_error_eV=analytic_err,commutator_max_error_eV=comm,residual_max_eV=resid,orthogonality_max_error=ortho,levels=levels,warnings=sorted(set(allwarnings)))
  reports.append(report)
  cases.append({'case':name,'group':group,'geometry_angstrom':[{'element':el,'position':pos.tolist()} for el,pos in atoms],'basis':[{'atom_index':a,'type':'s' if v is None else 'p','direction':None if v is None else v.tolist()} for a,v in basis],'edges':edges,'electrons':ne,'levels':refs[name]})
  return e,c,H,levels
 def star(el,group,preset,name,dirs):
  R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];dirs=np.array(dirs,float);atoms=[(el,np.zeros(3))]+[('H',R*u) for u in dirs];basis=[(0,None)]+[(0,u) for u in np.eye(3)]+[(i,None) for i in range(1,len(atoms))]
  es=P.onsite(el,'s');ep=P.onsite(el,'p');eh=P.onsite('H','s');ss=P.get(el,'H','ss_sigma',R,preset);sp=P.get(el,'H','sp_sigma',R,preset)
  if el=='O':
   g=np.sqrt(2)*P.element('O','s','H','s',dirs[0],R,preset);exp=np.linalg.eigvalsh([[es,0,g],[0,ep,-np.sqrt(2)*sp*dirs[0,2]],[g,-np.sqrt(2)*sp*dirs[0,2],eh]]).tolist()+E2(ep,eh,-np.sqrt(2)*sp*dirs[0,1])+[ep]
  elif el=='C':exp=E2(es,eh,2*ss)+E2(ep,eh,-2*sp/np.sqrt(3))*3
  else:exp=E2(es,eh,np.sqrt(3)*ss)+E2(ep,eh,-np.sqrt(1.5)*sp)*2+[ep]
  return test(name,group,atoms,basis,preset,6 if el=='B' else 8,exp)
 for preset,name in [('H2O_model','H2O model'),('H2O_lecture','H2O lecture')]:
  th=np.deg2rad(P.presets[preset]['geometry_overrides']['angle_degrees']);wd=[[0,np.sin(th/2),np.cos(th/2)],[0,-np.sin(th/2),np.cos(th/2)]];res=star('O','C2v',preset,name,wd)
  if 'lecture' in preset:assert abs(res[0][0]-P.onsite('O','s'))<1e-10 and abs(res[1][0,0])>1-1e-10
 cd=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3);bd=np.array([[np.cos(a),np.sin(a),0] for a in np.arange(3)*2*np.pi/3])
 ch=star('C','Td','CH4_model','CH4',cd);assert ch[3][0]['irrep']=='A1' and ch[3][1]['irrep']=='T2'
 bh=star('B','D3h','BH3_model','BH3',bd);assert [x for x in bh[3] if not x['electrons']][0]['irrep']=="A2''"
 axis=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
 for preset,name in [('SF6_model','SF6 model'),('SF6_lecture','SF6 lecture')]:
  R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];pp=P.get('F','F','pp_pi',np.sqrt(2)*R,preset);ps=P.get('F','F','pp_sigma',np.sqrt(2)*R,preset);tc=(pp-ps)/2;tt=-P.get('F','F','pp_sigma',2*R,preset);ef=P.onsite('F','p');la=ef+4*tc+tt;lt=ef-tt;le=ef-2*tc+tt
  exp=E2(P.onsite('S','s'),la,-np.sqrt(6)*P.get('S','F','sp_sigma',R,preset))+E2(P.onsite('S','p'),lt,-np.sqrt(2)*P.get('S','F','pp_sigma',R,preset))*3+[le]*2
  sf=test(name,'Oh',[('S',np.zeros(3))]+[('F',R*u) for u in axis],[(0,None)]+[(0,u) for u in np.eye(3)]+[(i,-axis[i-1]) for i in range(1,7)],preset,12,exp)
  assert sf[3][0]['irrep']=='A1g' and [x for x in sf[3] if x['electrons']][-1]['irrep']=='Eg'
  stage=np.array([[ef+tt,2*tc],[2*tc,ef+tt]]);assert np.ptp(np.linalg.eigvalsh(stage))>1e-8
 eh=P.onsite('H','s');preset='H2_fixture';R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];t=P.get('H','H','ss_sigma',R,preset)
 test('H2','D2h',[('H',np.array([0,0,-R/2])),('H',np.array([0,0,R/2]))],[(0,None),(1,None)],preset,2,[eh+t,eh-t])
 preset='H3_fixture';R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];t=P.get('H','H','ss_sigma',R,preset)
 for name,ne in [('H3+',2),('neutral H3',3)]:
  h3=test(name,'D3h',[('H',u*R/np.sqrt(3)) for u in bd],[(i,None) for i in range(3)],preset,ne,[eh+2*t]+[eh-t]*2)
  if ne==3:assert 'occupation_not_unique' in reports[-1]['warnings'] and 'Jahn_Teller_warning' in reports[-1]['warnings']
 preset='allyl_pi';R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];th=np.deg2rad(P.presets[preset]['geometry_overrides']['angle_degrees']);t=P.get('C','C','pp_pi',R,preset);ec=P.onsite('C','p')
 allyl=test('allyl cation pi','C2v',[('C',np.array([0,R*np.sin(th/2),R*np.cos(th/2)])),('C',np.zeros(3)),('C',np.array([0,-R*np.sin(th/2),R*np.cos(th/2)]))],[(i,np.array([1.,0,0])) for i in range(3)],preset,2,[ec+np.sqrt(2)*t,ec,ec-np.sqrt(2)*t])
 for preset,name,dirs,group,ne in [('sixH_lecture','six H Oh',axis,'Oh',6),('fourH_diagnostic','four H Td',cd,'Td',4)]:
  R=P.presets[preset]['geometry_overrides']['radius_angstrom']
  if group=='Oh':tcH=P.get('H','H','ss_sigma',np.sqrt(2)*R,preset);ttH=P.get('H','H','ss_sigma',2*R,preset);exp=[eh+4*tcH+ttH]+[eh-ttH]*3+[eh-2*tcH+ttH]*2
  else:t=P.get('H','H','ss_sigma',np.sqrt(8/3)*R,preset);exp=[eh+3*t]+[eh-t]*3
  test(name,group,[('H',R*u) for u in dirs],[(i,None) for i in range(len(dirs))],preset,ne,exp);assert 'occupation_not_unique' in reports[-1]['warnings']
 # Missing parameter lookups must error, not become zeros.
 try:P.get('N','N','pp_pi',1.,'H2_fixture')
 except ValueError:pass
 else:raise AssertionError('missing pair silently accepted')
 if check_reference:
  baseline=json.loads((ROOT/'evidence'/'reference_spectra.json').read_text());assert set(baseline)==set(refs)
  for name,levels in refs.items():
   assert len(levels)==len(baseline[name]),name
   for r,b in zip(levels,baseline[name]):
    assert r['irrep']==b['irrep'] and r['multiplicity']==b['multiplicity'] and r['electrons']==b['electrons'],name
    assert abs(r['energy_eV']-b['energy_eV'])<1e-9,(name,'committed reference changed',r,b)
 write_json(out/'checks'/'release1_validation.json',{'fixtures':reports,'SF6_combine_block_eV':stage.tolist(),'SF6_combine_split_eV':float(np.ptp(np.linalg.eigvalsh(stage))),'status':'PASS_CURRENT_FIXTURES_ONLY','reference_checked':check_reference,'warning_policy':'degenerate partially filled shell triggers occupation/JT/spin warning, not an assertion of an actual distorted ground state','missing_diatomics':['B2','C2','O2','F2','CO']})
 write_json(out/'checks'/'crystod_cases.json',cases)
 out.joinpath('checks').mkdir(exist_ok=True,parents=True)
 with (out/'checks'/'release1_energies.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 return refs
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=ROOT);p.add_argument('--groups',type=Path,default=ROOT/'groups');p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.data,a.out,a.groups)
