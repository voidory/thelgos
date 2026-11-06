import argparse,csv,json
from pathlib import Path
import numpy as np
from common import ROOT,write_json,csv_read
from parameters import Parameters

def E2(a,b,g):return np.linalg.eigvalsh([[a,g],[g,b]]).tolist()
def warn_shell(d,n):
 partial=0<n<2*d
 return ['occupation_not_unique','Jahn_Teller_warning','spin_state_not_determined'] if d>1 and partial else []
# [K1] Pair classification uses normalised positions, independent of bond length/radius.
def cistrans(p,q):
 return 'trans' if float(p@q)/(np.linalg.norm(p)*np.linalg.norm(q))<-.99 else 'cis'
# [K2] Irrep labels asserted against CONVENTIONS Appendix B, not only against the self-generated reference_spectra.json.
def expect_labels(name,levels,want):
 got=[x['irrep'] for x in levels]
 if got!=want:raise AssertionError(f'{name}: irrep sequence {got} != Appendix B {want}')
# [K14] Linear molecules run in a finite subgroup, where a pi pair spans two irreps (B2u+B3u in D2h, B1+B2 in C2v).
# Degenerate sets are then labelled by their exact irrep decomposition, and mapped to Dinfh/Cinfv names.
LINEAR_LABELS={'D2h':{'Ag':'sigma_g','B1u':'sigma_u','B2u+B3u':'pi_u','B2g+B3g':'pi_g'},'C2v':{'A1':'sigma','B1+B2':'pi'}}
def decompose(V,projectors,name):
 n={l:float(np.trace(V.T@pr@V)) for l,pr in projectors.items()};parts=sorted(l for l,x in n.items() if x>.5)
 if any(abs(n[l]-round(n[l]))>1e-3 for l in n) or abs(sum(n.values())-V.shape[1])>1e-6:raise AssertionError((name,'non-integer irrep decomposition',n))
 return '+'.join(parts)
# [K4] Level descriptor shared with the CrystOD comparison (see compare_crystod.py and CHECKER-CHANGES.md).
DESCRIPTOR_RULE=('character: sign of the bond energy sum 2*sum c_i c_j H_ij over central-ligand pairs (central topology) '
 'or all inter-atomic pairs (otherwise), averaged over the degenerate set; |x|<=1e-6 eV is nonbonding. '
 'dominant_fragment: central if the central-atom weight averaged over the degenerate set exceeds 0.5, else ligand; ligand_only when there is no central atom.')
def descriptors(H,c,basis,level_ix,topology,tol=1e-6):
 N=len(basis);central=topology=='central'
 scope=[(i,j) for i in range(N) for j in range(i+1,N) if basis[i][0]!=basis[j][0] and (not central or (basis[i][0]==0)!=(basis[j][0]==0))]
 out=[]
 for ix in level_ix:
  V=c[:,ix];d=len(ix);bond=sum(2*float(V[i]@V[j])*H[i,j] for i,j in scope)/d
  ch='bonding' if bond<-tol else ('antibonding' if bond>tol else 'nonbonding')
  if central:frag='central' if sum(float(V[i]@V[i]) for i in range(N) if basis[i][0]==0)/d>.5 else 'ligand'
  else:frag='ligand_only'
  out.append({'character':ch,'dominant_fragment':frag,'bond_energy_eV':float(bond)})
 return out

def run(data,out,groups,check_reference=True):
 P=Parameters(data);reports=[];rows=[];refs={};cases=[]
 def test(name,group,atoms,basis,preset,ne,expected,topology='all',linear=False):
  p=P.presets[preset];N=len(basis);H=np.zeros((N,N));edges=[]
  for i,(a,v) in enumerate(basis):H[i,i]=P.onsite(atoms[a][0],'s' if v is None else 'p')
  for a in range(len(atoms)):
   for b in range(a+1,len(atoms)):
    u=atoms[b][1]-atoms[a][1];d=float(np.linalg.norm(u));u/=d
    # [K1] Explicit topology per fixture replaces matching on case names, preset ids and central-element lists.
    if topology=='chain':typ='nearest_neighbor' if b-a==1 else 'absent'
    elif topology=='central':typ='central_ligand' if a==0 else cistrans(atoms[a][1],atoms[b][1])
    elif topology=='cistrans':typ=cistrans(atoms[a][1],atoms[b][1])
    elif topology=='all':typ='all'
    else:raise ValueError(f'{name}: unknown topology {topology!r}')
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
  k=0;left=ne;levels=[];allwarnings=[];level_ix=[]
  while k<N:
   ix=np.flatnonzero(abs(e-e[k])<1e-8);ix=ix[ix>=k];level_ix.append(ix);V=c[:,ix];weights={l:float(np.trace(V.T@pr@V)/len(ix)) for l,pr in projectors.items()};l=max(weights,key=weights.get)
   if linear and weights[l]<=.999:l=decompose(V,projectors,name)
   assert weights.get(l,1.)>.999 or '+' in l,(name,weights)
   d=len(ix);n=min(left,2*d);left-=n;warnings=warn_shell(d,n);allwarnings.extend(warnings)
   row=dict(case=name,irrep=l,energy_eV=float(e[k]),multiplicity=d,electrons=n,warnings=';'.join(warnings));rows.append(row);levels.append(row);k=int(ix[-1])+1
  assert left==0
  resid=float(np.linalg.norm(H@c-c*e,axis=0).max());ortho=float(np.max(abs(c.T@c-np.eye(N))));assert resid<1e-10 and ortho<1e-10
  refs[name]=[dict(irrep=x['irrep'],energy_eV=x['energy_eV'],multiplicity=x['multiplicity'],electrons=x['electrons']) for x in levels]
  report=dict(case=name,status='PASS',analytic_max_error_eV=analytic_err,commutator_max_error_eV=comm,residual_max_eV=resid,orthogonality_max_error=ortho,levels=levels,warnings=sorted(set(allwarnings)))
  reports.append(report)
  desc=descriptors(H,c,basis,level_ix,topology)
  cases.append({'case':name,'group':group,'topology':topology,'descriptor_rule':DESCRIPTOR_RULE,'geometry_angstrom':[{'element':el,'position':pos.tolist()} for el,pos in atoms],'basis':[{'atom_index':a,'type':'s' if v is None else 'p','direction':None if v is None else v.tolist()} for a,v in basis],'edges':edges,'electrons':ne,'levels':[dict(r,descriptor=dd) for r,dd in zip(refs[name],desc)]})
  return e,c,H,levels
 def star(el,group,preset,name,dirs):
  R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];dirs=np.array(dirs,float);atoms=[(el,np.zeros(3))]+[('H',R*u) for u in dirs];basis=[(0,None)]+[(0,u) for u in np.eye(3)]+[(i,None) for i in range(1,len(atoms))]
  es=P.onsite(el,'s');ep=P.onsite(el,'p');eh=P.onsite('H','s');ss=P.get(el,'H','ss_sigma',R,preset);sp=P.get(el,'H','sp_sigma',R,preset)
  if el=='O':
   g=np.sqrt(2)*P.element('O','s','H','s',dirs[0],R,preset);exp=np.linalg.eigvalsh([[es,0,g],[0,ep,-np.sqrt(2)*sp*dirs[0,2]],[g,-np.sqrt(2)*sp*dirs[0,2],eh]]).tolist()+E2(ep,eh,-np.sqrt(2)*sp*dirs[0,1])+[ep]
  elif el=='C':exp=E2(es,eh,2*ss)+E2(ep,eh,-2*sp/np.sqrt(3))*3
  else:exp=E2(es,eh,np.sqrt(3)*ss)+E2(ep,eh,-np.sqrt(1.5)*sp)*2+[ep]
  return test(name,group,atoms,basis,preset,6 if el=='B' else 8,exp,topology='central')
 for preset,name in [('H2O_model','H2O model'),('H2O_lecture','H2O lecture')]:
  th=np.deg2rad(P.presets[preset]['geometry_overrides']['angle_degrees']);wd=[[0,np.sin(th/2),np.cos(th/2)],[0,-np.sin(th/2),np.cos(th/2)]];res=star('O','C2v',preset,name,wd)
  if 'lecture' in preset:assert abs(res[0][0]-P.onsite('O','s'))<1e-10 and abs(res[1][0,0])>1-1e-10
  expect_labels(name,res[3],['A1','B2','A1','B1','A1','B2'])
  assert [x['irrep'] for x in res[3] if abs(x['energy_eV']-P.onsite('O','p'))<1e-10]==['B1'],(name,'out-of-plane O px must be the B1 nonbonding level')
 cd=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3);bd=np.array([[np.cos(a),np.sin(a),0] for a in np.arange(3)*2*np.pi/3])
 ch=star('C','Td','CH4_model','CH4',cd);assert ch[3][0]['irrep']=='A1' and ch[3][1]['irrep']=='T2';expect_labels('CH4',ch[3],['A1','T2','A1','T2'])
 bh=star('B','D3h','BH3_model','BH3',bd);assert [x for x in bh[3] if not x['electrons']][0]['irrep']=="A2''";expect_labels('BH3',bh[3],["A1'","E'","A2''","A1'","E'"])
 axis=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
 for preset,name in [('SF6_model','SF6 model'),('SF6_lecture','SF6 lecture')]:
  R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];pp=P.get('F','F','pp_pi',np.sqrt(2)*R,preset);ps=P.get('F','F','pp_sigma',np.sqrt(2)*R,preset);tc=(pp-ps)/2;tt=-P.get('F','F','pp_sigma',2*R,preset);ef=P.onsite('F','p');la=ef+4*tc+tt;lt=ef-tt;le=ef-2*tc+tt
  exp=E2(P.onsite('S','s'),la,-np.sqrt(6)*P.get('S','F','sp_sigma',R,preset))+E2(P.onsite('S','p'),lt,-np.sqrt(2)*P.get('S','F','pp_sigma',R,preset))*3+[le]*2
  sf=test(name,'Oh',[('S',np.zeros(3))]+[('F',R*u) for u in axis],[(0,None)]+[(0,u) for u in np.eye(3)]+[(i,-axis[i-1]) for i in range(1,7)],preset,12,exp,topology='central')
  assert sf[3][0]['irrep']=='A1g' and [x for x in sf[3] if x['electrons']][-1]['irrep']=='Eg';expect_labels(name,sf[3],['A1g','T1u','Eg','A1g','T1u'])
  stage=np.array([[ef+tt,2*tc],[2*tc,ef+tt]]);assert np.ptp(np.linalg.eigvalsh(stage))>1e-8
 eh=P.onsite('H','s');preset='H2_fixture';R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];t=P.get('H','H','ss_sigma',R,preset)
 h2=test('H2','D2h',[('H',np.array([0,0,-R/2])),('H',np.array([0,0,R/2]))],[(0,None),(1,None)],preset,2,[eh+t,eh-t]);expect_labels('H2',h2[3],['Ag','B1u'])
 preset='H3_fixture';R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];t=P.get('H','H','ss_sigma',R,preset)
 for name,ne in [('H3+',2),('neutral H3',3)]:
  h3=test(name,'D3h',[('H',u*R/np.sqrt(3)) for u in bd],[(i,None) for i in range(3)],preset,ne,[eh+2*t]+[eh-t]*2)
  expect_labels(name,h3[3],["A1'","E'"])
  if ne==3:assert 'occupation_not_unique' in reports[-1]['warnings'] and 'Jahn_Teller_warning' in reports[-1]['warnings']
 preset='allyl_pi';R=P.presets[preset]['geometry_overrides']['bond_length_angstrom'];th=np.deg2rad(P.presets[preset]['geometry_overrides']['angle_degrees']);t=P.get('C','C','pp_pi',R,preset);ec=P.onsite('C','p')
 allyl=test('allyl cation pi','C2v',[('C',np.array([0,R*np.sin(th/2),R*np.cos(th/2)])),('C',np.zeros(3)),('C',np.array([0,-R*np.sin(th/2),R*np.cos(th/2)]))],[(i,np.array([1.,0,0])) for i in range(3)],preset,2,[ec+np.sqrt(2)*t,ec,ec-np.sqrt(2)*t],topology='chain')
 expect_labels('allyl cation pi',allyl[3],['B1','A2','B1'])
 occ=[x for x in allyl[3] if x['electrons']];vac=[x for x in allyl[3] if not x['electrons']]
 assert occ[-1]['irrep']=='B1' and vac[0]['irrep']=='A2',('allyl cation: HOMO must be B1, LUMO A2 (spec F2.10)',occ[-1],vac[0])
 for preset,name,dirs,group,ne in [('sixH_lecture','six H Oh',axis,'Oh',6),('fourH_diagnostic','four H Td',cd,'Td',4)]:
  R=P.presets[preset]['geometry_overrides']['radius_angstrom']
  if group=='Oh':tcH=P.get('H','H','ss_sigma',np.sqrt(2)*R,preset);ttH=P.get('H','H','ss_sigma',2*R,preset);exp=[eh+4*tcH+ttH]+[eh-ttH]*3+[eh-2*tcH+ttH]*2
  else:t=P.get('H','H','ss_sigma',np.sqrt(8/3)*R,preset);exp=[eh+3*t]+[eh-t]*3
  lig=test(name,group,[('H',R*u) for u in dirs],[(i,None) for i in range(len(dirs))],preset,ne,exp,topology='cistrans' if group=='Oh' else 'all');assert 'occupation_not_unique' in reports[-1]['warnings']
  expect_labels(name,lig[3],['A1g','T1u','Eg'] if group=='Oh' else ['A1','T2'])
 # [K14] Release 1 diatomics (gate C2). Geometry: NIST CCCBDB r_e (evidence/cccbdb_*.html); atom 1 at -R/2, atom 2 at +R/2 on z.
 # Explicit analytic blocks, independent of P.block: homonuclear D2h sigma_g D(es+Vss, ep-Vpps; -Vsp), sigma_u D(es-Vss, ep+Vpps; Vsp),
 # pi_u ep+Vppp (x2), pi_g ep-Vppp (x2). CO (C2v): 4x4 sigma block [sC,zC,sO,zO] and pi D(epC, epO; Vppp) (x2).
 # Reference orderings (CONVENTIONS A.6): B2, C2 pi_u below sigma_g(2p); O2, F2 sigma_g(2p) below pi_u. A disagreement is recorded as an
 # Orbital-ordering fidelity warning, never tuned away.
 diatomic_reports=[]
 REF_ORDER={'B2':'pi_u_below_sigma_g_2p','C2':'pi_u_below_sigma_g_2p','N2':'pi_u_below_sigma_g_2p','O2':'sigma_g_2p_below_pi_u','F2':'sigma_g_2p_below_pi_u'}
 # Ground states as listed by NIST CCCBDB (same evidence pages): an open shell is expected exactly when the reference is not a singlet.
 REF_STATE={'B2':('3Sigma_g',3),'C2':('1Sigma_g+',1),'N2':('1Sigma_g',1),'O2':('3Sigma_g',3),'F2':('1Sigma_g',1),'CO':('1Sigma',1)}
 for name,a,b,ne in [('B2','B','B',6),('C2','C','C',8),('N2','N','N',10),('O2','O','O',12),('F2','F','F',14),('CO','C','O',10)]:
  preset=f'{name}_model';R=P.presets[preset]['geometry_overrides']['bond_length_angstrom']
  atoms=[(a,np.array([0,0,-R/2])),(b,np.array([0,0,R/2]))];basis=[(i,v) for i in range(2) for v in [None]+list(np.eye(3))]
  Vss,Vsp,Vps,Vpp=[P.get(a,b,ch,R,preset) for ch in ['ss_sigma','sp_sigma','pp_sigma','pp_pi']]
  if a==b:
   es,ep=P.onsite(a,'s'),P.onsite(a,'p');group='D2h'
   exp=E2(es+Vss,ep-Vps,-Vsp)+E2(es-Vss,ep+Vps,Vsp)+[ep+Vpp]*2+[ep-Vpp]*2
  else:
   group='C2v';ea,pa,eb,pb=P.onsite(a,'s'),P.onsite(a,'p'),P.onsite(b,'s'),P.onsite(b,'p')
   sig=np.array([[ea,0,Vss,Vsp],[0,pa,-Vsp,Vps],[Vss,-Vsp,eb,0],[Vsp,Vps,0,pb]])
   exp=np.linalg.eigvalsh(sig).tolist()+E2(pa,pb,Vpp)*2
  e,c,H,lv=test(name,group,atoms,basis,preset,ne,exp,topology='all',linear=True)
  lin=[LINEAR_LABELS[group][x['irrep']] for x in lv]
  if group=='D2h':assert sorted(set(x['irrep'] for x in lv))==['Ag','B1u','B2g+B3g','B2u+B3u'],(name,lv)
  else:assert sorted(set(x['irrep'] for x in lv))==['A1','B1+B2'],(name,lv)
  desc=cases[-1]['levels'];occ=[k for k,x in enumerate(lv) if x['electrons']]
  bonding=sum(x['electrons'] for x,d in zip(lv,desc) if d['descriptor']['character']=='bonding');anti=sum(x['electrons'] for x,d in zip(lv,desc) if d['descriptor']['character']=='antibonding')
  rep={'case':name,'group':group,'bond_length_angstrom':R,'linear_labels':lin,'counting_bond_order':(bonding-anti)/2,'fidelity_warnings':[]}
  if name in REF_ORDER:
   sg=[x['energy_eV'] for x,l in zip(lv,lin) if l=='sigma_g'];pu=[x['energy_eV'] for x,l in zip(lv,lin) if l=='pi_u'][0]
   model='pi_u_below_sigma_g_2p' if pu<max(sg) else 'sigma_g_2p_below_pi_u'
   rep['orbital_ordering']={'reference':REF_ORDER[name],'model':model,'sigma_g_2p_taken_as':'higher of the two sigma_g levels','agrees':model==REF_ORDER[name]}
   if model!=REF_ORDER[name]:rep['fidelity_warnings'].append('orbital_ordering_disagrees_with_reference')
  st,mult=REF_STATE[name];open_shell=bool(reports[-1]['warnings'])
  rep['spin_state']={'reference_ground_state':st,'model_open_shell':open_shell,'agrees':open_shell==(mult>1)}
  if open_shell!=(mult>1):rep['fidelity_warnings'].append('spin_state_disagrees_with_reference')
  if name=='O2':assert 'occupation_not_unique' in reports[-1]['warnings'] and lin[occ[-1]]=='pi_g',('O2 must end in a partly filled pi_g',lv)
  if name=='CO':
   h=occ[-1];v=c[:,[i for i in range(len(e)) if abs(e[i]-lv[h]['energy_eV'])<1e-8]];wC=float((v[:4]**2).sum()/v.shape[1])
   assert lin[h]=='sigma' and wC>.5 and rep['counting_bond_order']==3,('CO HOMO / bond order',lin[h],wC,rep['counting_bond_order'])
   rep['HOMO_carbon_weight']=wC
  diatomic_reports.append(rep)
 # Missing parameter lookups must error, not become zeros.
 # [K16] v1.2 used N-N as the absent element; N is now parameterised (decision B4), so neon stands in.
 try:P.get('Ne','Ne','pp_pi',1.,'H2_fixture')
 except ValueError:pass
 else:raise AssertionError('missing pair silently accepted')
 # [K3, revised in patch 3] An existing element pair with an absent channel: O-H has ss/sp rows but no pp_pi row.
 assert any(sorted([r['element_a'],r['element_b']])==['H','O'] for r in P.rows) and not any(sorted([r['element_a'],r['element_b']])==['H','O'] and r['channel']=='pp_pi' for r in P.rows),'K3 precondition: O-H rows exist without pp_pi'
 try:P.get('O','H','pp_pi',1.,'H2O_model')
 except (ValueError,LookupError):pass
 else:raise AssertionError('absent O-H pp_pi channel silently accepted')
 if check_reference:
  baseline=json.loads((ROOT/'evidence'/'reference_spectra.json').read_text());assert set(baseline)==set(refs)
  for name,levels in refs.items():
   assert len(levels)==len(baseline[name]),name
   for r,b in zip(levels,baseline[name]):
    assert r['irrep']==b['irrep'] and r['multiplicity']==b['multiplicity'] and r['electrons']==b['electrons'],name
    assert abs(r['energy_eV']-b['energy_eV'])<1e-9,(name,'committed reference changed',r,b)
 write_json(out/'checks'/'release1_validation.json',{'fixtures':reports,'SF6_combine_block_eV':stage.tolist(),'SF6_combine_split_eV':float(np.ptp(np.linalg.eigvalsh(stage))),'status':'PASS_CURRENT_FIXTURES_ONLY','reference_checked':check_reference,'warning_policy':'degenerate partially filled shell triggers occupation/JT/spin warning, not an assertion of an actual distorted ground state','diatomics':diatomic_reports,'missing_diatomics':[]})
 write_json(out/'checks'/'crystod_cases.json',cases)
 out.joinpath('checks').mkdir(exist_ok=True,parents=True)
 with (out/'checks'/'release1_energies.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 return refs
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=ROOT);p.add_argument('--groups',type=Path,default=ROOT/'groups');p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.data,a.out,a.groups)
