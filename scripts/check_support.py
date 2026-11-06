import csv,json,tempfile,shutil
from pathlib import Path
import numpy as np
from markdown_it import MarkdownIt
from jsonschema import Draft202012Validator,ValidationError
from pysktb._params import get_hop_int
from common import ROOT,write_json,sha
from parameters import Parameters
from check_release1 import run as run_physics

def run(out):
 P=Parameters(ROOT)
 # Untouched data files must match committed independent-reference spectra; a pair mutation must not.
 with tempfile.TemporaryDirectory() as tmp:
  root=Path(tmp);shutil.copytree(ROOT/'params',root/'params');rows=list(csv.DictReader((root/'params/pairs.csv').open()))
  row=next(r for r in rows if sorted([r['element_a'],r['element_b']])==['H','O'] and r['channel']=='sp_sigma');old=row['eta'];row['eta']=str(float(old)*1.01)
  with (root/'params/pairs.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  try:run_physics(root,root/'out',ROOT/'groups')
  except AssertionError as err:detected='committed reference changed' in str(err);reason='reference spectrum mismatch' if detected else str(err)
  else:detected=False;reason='not detected'
  assert detected,reason
 write_json(out/'checks/parameter_mutation.json',{'status':'PASS','mutation':'O-H sp_sigma eta multiplied by 1.01 in temporary pairs.csv copy','unchanged_committed_data':True,'detected':detected,'mechanism':reason})
 data=json.loads((ROOT/'params/presets.json').read_text());schema=json.loads((ROOT/'params/presets.schema.json').read_text());Draft202012Validator(schema).validate(data)
 bad=json.loads(json.dumps(data));del bad['presets'][0]['edge_policy']
 try:Draft202012Validator(schema).validate(bad)
 except ValidationError:pass
 else:raise AssertionError('schema failed to reject missing edge policy')
 write_json(out/'checks/preset_schema_validation.json',{'status':'PASS','valid_presets':len(data['presets']),'missing_edge_policy_rejected':True})
 # Unscaled values come from actual pair rows, cancelling their declared single H scale.
 specs=[('O-H sp at water','O','H','sp_sigma',.958,'H2O_model'),('O-H ss at water','O','H','ss_sigma',.958,'H2O_model'),('C-H ss at methane','C','H','ss_sigma',1.087,'CH4_model'),('C-H sp at methane','C','H','sp_sigma',1.087,'CH4_model'),('B-H ss at borane','B','H','ss_sigma',1.19,'BH3_model'),('B-H sp at borane','B','H','sp_sigma',1.19,'BH3_model'),('H2','H','H','ss_sigma',1.,'H2_fixture'),('H3','H','H','ss_sigma',1.,'H3_fixture'),('six H cis','H','H','ss_sigma',np.sqrt(2),'sixH_lecture'),('six H trans','H','H','ss_sigma',2.,'sixH_lecture'),('four H edge','H','H','ss_sigma',1.087*np.sqrt(8/3),'fourH_diagnostic'),('allyl pi','C','C','pp_pi',1.4,'allyl_pi')]
 rows=[]
 for case,a,b,ch,d,preset in specs:
  effective=P.get(a,b,ch,d,preset);scale=P.presets[preset]['hydrogen_scale'] if 'H' in [a,b] else 1.;rows.append(dict(case=case,channel=ch,d_angstrom=d,unscaled_eV=effective/scale,effective_eV=effective,scale=scale))
 assert abs(rows[0]['unscaled_eV']-15.277130068296428)<1e-12
 out.joinpath('checks').mkdir(parents=True,exist_ok=True)
 with (out/'checks/unscaled_reference_couplings.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 # Independent angular implementation: full polynomial/tensor contractions, not copies of the SK table.
 r3=np.sqrt(3);T=np.array([np.diag([-.5,-.5,1]),[[0,0,r3/2],[0,0,0],[r3/2,0,0]],[[0,0,0],[0,0,r3/2],[0,r3/2,0]],np.diag([r3/2,-r3/2,0]),[[0,r3/2,0],[r3/2,0,0],[0,0,0]]]);order=[0,1,2,3,8,6,5,7,4];params=dict(V_sss=-1.4,V_sps=1.84,V_pps=3.24,V_ppp=-.81,V_sds=-3.16,V_pds=-2.95,V_pdp=1.36)
 maxerr=0.;reversal=0.;count=0
 for a in [-2,-1,0,1,2]:
  for b in [-2,-1,0,1,2]:
   for c in [-2,-1,0,1,2]:
    if (a,b,c)==(0,0,0):continue
    u=np.array([a,b,c],float);u/=np.linalg.norm(u);raw=get_hop_int(**params,l=u[0],m=u[1],n=u[2]);S=np.array([[raw[i][j] or 0 for j in order] for i in order]);f=np.einsum('i,kij,j->k',u,T,u)
    A=np.zeros((9,9));A[0,0]=params['V_sss'];A[0,1:4]=u*params['V_sps'];A[1:4,0]=-u*params['V_sps'];A[1:4,1:4]=params['V_ppp']*np.eye(3)+(params['V_pps']-params['V_ppp'])*np.outer(u,u);A[0,4:]=f*params['V_sds'];A[4:,0]=A[0,4:];A[1:4,4:]=np.outer(u,f)*params['V_pds']+(2/r3)*(np.einsum('kij,j->ik',T,u)-np.outer(u,f))*params['V_pdp'];A[4:,1:4]=-A[1:4,4:].T
    maxerr=max(maxerr,float(np.max(abs(S-A))));rawB=get_hop_int(**params,l=-u[0],m=-u[1],n=-u[2]);B=np.array([[rawB[i][j] or 0 for j in order] for i in order]);reversal=max(reversal,float(np.max(abs(B-S.T))));assert abs(np.sum(f*f)-1)<1e-12;count+=1
 assert maxerr<1e-12 and reversal<1e-12
 write_json(out/'checks/angular_validation.json',{'status':'PASS','reference':'pysktb 0.5.6 get_hop_int','directions':count,'max_error':maxerr,'reverse_transpose_error':reversal,'channels':'s/p, sd, pd; dd not checked'})
 tokens=MarkdownIt('commonmark').enable('table').parse((ROOT/'CONVENTIONS.md').read_text());tables=[];cur=None
 for t in tokens:
  if t.type=='table_open':cur=[]
  elif t.type=='tr_open' and cur is not None:cur.append(0)
  elif t.type in ['th_open','td_open'] and cur is not None:cur[-1]+=1
  elif t.type=='table_close':tables.append(cur);cur=None
 assert tables and tables[0]==[2]*11,tables[:1];assert all(len(set(t))==1 for t in tables)
 write_json(out/'checks/markdown_validation.json',{'status':'PASS','sign_table_columns_by_row':tables[0],'rectangular_tables':len(tables)})
 lock=(ROOT/'requirements-checker.lock').read_text();blocks=[b for b in lock.split('\n\n') if '==' in b];assert all('--hash=sha256:' in b for b in blocks)
 write_json(out/'checks/lock_validation.json',{'status':'PASS','hashed_distributions':len(blocks),'scope':'complete CPython3.12/Linux checker dependency graph; bootstrap pip excluded','clean_hash_install':'performed in separate environment; see evidence/hash_install.txt','production_python_pin':'3.12.7; this CPython check is not the browser gate','runtime_python':'printed at run time, not committed (CHECKER-CHANGES.md K12)'})
 # [K12] The host interpreter version is reported, not written into a byte-compared committed file.
 print('check_support: runtime Python',__import__('platform').python_version())
 # Checker comparator unit tests are mock protocol tests only, not claimed CrystOD results.
 from compare_crystod import compare,run as crystod_run
 import copy
 # [K10] Mock inputs carry K4 descriptors and K7 contract-basis weights; mutations are split into blocking and order-review classes (K11).
 Dl=lambda ch,fr:{'character':ch,'dominant_fragment':fr}
 case={'case':'protocol_unit_test','group':'C2v','geometry_angstrom':[],'basis':[],'levels':[{'irrep':'A1','energy_eV':-2.,'multiplicity':1,'descriptor':Dl('bonding','ligand')},{'irrep':'A1','energy_eV':-1.,'multiplicity':1,'descriptor':Dl('antibonding','central')}]}
 sample={'case':'protocol_unit_test','route':'extended_huckel','pyscf_used':False,'axis_map_checker_to_contract':np.eye(3).tolist(),'basis_description':'mock','basis_differences':'mock protocol test only','irrep_map_to_contract':{'a1':'A1'},'point_group_full':'C2v','level_identity_evidence':'mock descriptors assigned from mock eigenvectors','levels':[{'irrep':'a1','energy':-5.,'multiplicity':1,'canonical_level_id':'extra','contract_basis_weight':.1,'descriptor':Dl('bonding','ligand')},{'irrep':'a1','energy':-2.,'multiplicity':1,'canonical_level_id':'A1:1','contract_basis_weight':.9,'descriptor':Dl('bonding','ligand')},{'irrep':'a1','energy':-1.,'multiplicity':1,'canonical_level_id':'A1:2','contract_basis_weight':.9,'descriptor':Dl('antibonding','central')}]}
 good=compare(case,sample);assert good['status']=='PASS' and len(good['excluded_outside_contract_basis'])==1,good
 for mode in ['point_group','degeneracy','labels','extra_inside_contract_basis','missing_weight','missing_descriptor']:
  bad=copy.deepcopy(sample)
  if mode=='point_group':bad['point_group_full']='Td'
  elif mode=='degeneracy':bad['levels'][1]['multiplicity']=2
  elif mode=='labels':bad['irrep_map_to_contract']={'a1':'B1'}
  elif mode=='extra_inside_contract_basis':bad['levels'][0]['contract_basis_weight']=.9
  elif mode=='missing_weight':del bad['levels'][1]['contract_basis_weight']
  else:del bad['levels'][1]['descriptor']
  try:rejected=compare(case,bad)['status']=='DISAGREEMENT_ORACLES_PREVAIL'
  except ValueError:rejected=True
  assert rejected,mode
 for mode in ['energy_swap','descriptor_swap']:
  bad=copy.deepcopy(sample)
  if mode=='energy_swap':bad['levels'][1]['energy']=0.
  else:bad['levels'][1]['descriptor'],bad['levels'][2]['descriptor']=bad['levels'][2]['descriptor'],bad['levels'][1]['descriptor']
  assert compare(case,bad)['status']=='ORDER_REVIEW_FOR_OWNER',mode
 with tempfile.TemporaryDirectory() as temp:
  p=Path(temp);write_json(p/'checks/crystod_manifest.json',{'status':'PASS'})
  assert crystod_run(p,cases_file=out/'checks/crystod_cases.json')==2
  assert json.loads((p/'checks/crystod_manifest.json').read_text())['status']=='BLOCKED'
 write_json(out/'checks/crystod_protocol_tests.json',{'status':'PASS_MOCK_PROTOCOL_ONLY','blocking_mutations_rejected':['point_group','degeneracy','labels','extra_inside_contract_basis','missing_weight','missing_descriptor'],'order_mutations_sent_to_review':['energy_swap','descriptor_swap'],'outside_contract_level_excluded':True,'hand_set_PASS_ignored':True,'real_archive_or_adapter_executed':False})
