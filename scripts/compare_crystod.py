"""Never trusts a committed PASS. Computes archive identity and comparison status each run.
The real archive/API is unavailable. The adapter protocol is specified, not asserted tested with CrystOD.
"""
import argparse,hashlib,json,re,subprocess,sys,tomllib,zipfile
from pathlib import Path
import numpy as np
from common import ROOT,write_json,sha

def version_from_zip(path):
 with zipfile.ZipFile(path) as z:
  for name in sorted(z.namelist()):
   if name.endswith('.dist-info/METADATA') or name.endswith('PKG-INFO'):
    raw=z.read(name).decode('utf-8',errors='replace');pkg=re.search(r'^Name:\s*(.+)$',raw,re.M);v=re.search(r'^Version:\s*(.+)$',raw,re.M)
    if pkg and 'crystod' in pkg.group(1).lower() and v:return v.group(1).strip(),name
   if name.endswith('pyproject.toml'):
    d=tomllib.loads(z.read(name).decode());project=d.get('project',d.get('tool',{}).get('poetry',{}))
    if 'crystod' in project.get('name','').lower() and project.get('version'):return str(project['version']),name
   if name.lower().endswith(('version.txt','_version.py','__init__.py')) and 'crystod' in name.lower():
    txt=z.read(name).decode(errors='replace');m=re.search(r"__version__\s*=\s*['\"]([^'\"]+)['\"]",txt)
    if m:return m.group(1),name
 raise ValueError('No unambiguous embedded CrystOD version: inspect archive and extend extractor with source-backed rule')

# [K5] Explicit failure messages; asserts vanish under python -O and gave empty "AssertionError: " reasons.
def require(cond,msg):
 if not cond:raise ValueError(msg)

# [K6] Full point group of linear molecules from geometry, not from the case name.
def full_point_group(case):
 els=[a['element'] for a in case['geometry_angstrom']];pos=np.array([a['position'] for a in case['geometry_angstrom']],float);cen=pos-pos.mean(0)
 if len(pos)>=2 and np.linalg.matrix_rank(cen,tol=1e-8)<=1:
  centro=all(any(e2==e1 and np.linalg.norm(p2+p1)<1e-8 for e2,p2 in zip(els,cen)) for e1,p1 in zip(els,cen))
  return 'Dinfh' if centro else 'Cinfv'
 return case['group']

# [K7] Levels the checker places mostly outside the contract basis are excluded before counting, not treated as disagreements.
CONTRACT_WEIGHT_MIN=0.5
# [K8] Within-irrep order is compared through eigenvector-derived descriptors, not adapter-chosen level ids.
DESCRIPTOR_KEYS=('character','dominant_fragment')
def dkey(level):return [level['descriptor'][k] for k in DESCRIPTOR_KEYS]

def compare(case,result):
 c=case['case']
 require(result.get('route')=='extended_huckel' and result.get('pyscf_used') is False,f'{c}: adapter must report route=extended_huckel and pyscf_used=false')
 require(result.get('case')==c,f"{c}: adapter answered for case {result.get('case')!r}")
 Q=np.array(result['axis_map_checker_to_contract'],float);require(Q.shape==(3,3) and np.max(abs(Q.T@Q-np.eye(3)))<1e-10,f'{c}: axis_map_checker_to_contract is not a 3x3 orthogonal matrix')
 require(result.get('basis_description') and result.get('basis_differences'),f'{c}: basis_description and basis_differences are mandatory')
 # Explicit label map must be backed by the adapter's actual coordinate convention.
 mapping=result.get('irrep_map_to_contract');require(isinstance(mapping,dict) and mapping,f'{c}: irrep_map_to_contract missing or empty')
 got=result['levels'];unmapped=sorted({x['irrep'] for x in got if x['irrep'] not in mapping});require(not unmapped,f'{c}: checker irreps without a contract mapping: {unmapped}')
 for x in got:
  w=x.get('contract_basis_weight');require(isinstance(w,(int,float)) and 0<=w<=1,f"{c}: level {x.get('canonical_level_id')!r} lacks contract_basis_weight in [0,1]")
  require(all(k in x.get('descriptor',{}) for k in DESCRIPTOR_KEYS),f"{c}: level {x.get('canonical_level_id')!r} lacks descriptor keys {DESCRIPTOR_KEYS}")
 expected=case['levels'];require(all('descriptor' in x for x in expected),f'{c}: crystod_cases.json predates descriptors; regenerate it with check_release1.py')
 inside=[x for x in got if x['contract_basis_weight']>=CONTRACT_WEIGHT_MIN];outside=[x for x in got if x['contract_basis_weight']<CONTRACT_WEIGHT_MIN]
 ours_irreps={x['irrep'] for x in expected};shared=set(mapping[x['irrep']] for x in inside)&ours_irreps;require(shared,f'{c}: no irrep shared between checker (inside contract basis) and contract')
 full_expected=full_point_group(case);point_ok=result['point_group_full']==full_expected
 checks=[]
 for label in sorted(shared):
  ours=sorted([x for x in expected if x['irrep']==label],key=lambda x:x['energy_eV']);theirs=sorted([x for x in inside if mapping[x['irrep']]==label],key=lambda x:x['energy'])
  od=[dkey(x) for x in ours];td=[dkey(x) for x in theirs]
  dims=[x['multiplicity'] for x in theirs];want=[x['multiplicity'] for x in ours]
  checks.append({'irrep':label,'shared_basis_counts_match':len(ours)==len(theirs),'degeneracy_pattern_matches':dims==want,'label_after_axis_mapping':label,'descriptor_order_matches':od==td,'expected_descriptor_order':od,'checker_descriptor_order':td,'checker_canonical_level_ids':[x.get('canonical_level_id') for x in theirs],'identity_assignment_evidence':result.get('level_identity_evidence')})
 ok=point_ok and all(k['shared_basis_counts_match'] and k['degeneracy_pattern_matches'] and k['descriptor_order_matches'] for k in checks)
 return {'case':c,'geometry_angstrom':case['geometry_angstrom'],'contract_basis':case['basis'],'checker_basis':result['basis_description'],'basis_differences':result['basis_differences'],'axis_map_checker_to_contract':Q.tolist(),'irrep_map':mapping,'point_group_matches':point_ok,'expected_full_point_group':full_expected,'checker_full_point_group':result['point_group_full'],'shared_irrep_checks':checks,'contract_irreps_not_shared':sorted(ours_irreps-shared),'excluded_outside_contract_basis':[{'irrep':x['irrep'],'canonical_level_id':x.get('canonical_level_id'),'contract_basis_weight':x['contract_basis_weight']} for x in outside],'status':'PASS' if ok else 'DISAGREEMENT_ORACLES_PREVAIL','energy_comparison':'NOT_PERFORMED; raw energies used only to sort within each irrep','composition_and_bond_orders':'NOT_COMPARED'}

def run(out,archive=None,adapter=None,cases_file=None):
 cases=json.loads((cases_file or ROOT/'checks/crystod_cases.json').read_text())
 report={'schema_version':'1.2','status':'BLOCKED','archive_filename':None,'archive_sha256':None,'version':None,'comparisons':[],'allowed':['point_group','shared_irrep_degeneracy','axis_mapped_labels','within_irrep_descriptor_order'],'comparison_rules':{'contract_basis_weight_min':CONTRACT_WEIGHT_MIN,'descriptor_keys':list(DESCRIPTOR_KEYS),'linear_point_group':'Dinfh/Cinfv from geometry','changes':'CHECKER-CHANGES.md K5-K9'},'prohibited':['energy_values','compositions','bond_orders'],'precedence':'Appendix A oracles win; discrepancies recorded for the owner','engine_UI_gate':'separate, not exercised by this script','diatomic_gate':'B2,C2,O2,F2,CO missing; N2 scope undecided','required_cases':[x['case'] for x in cases]+['B2','C2','O2','F2','CO']}
 if archive is None or not archive.is_file():report['reason']='Original supplied CrystOD zip unavailable; cannot compute version/hash or execute checker';code=2
 else:
  report.update(archive_filename=archive.name,archive_sha256=sha(archive))
  try:
   version,source=version_from_zip(archive);report.update(version=version,version_metadata_member=source)
   if adapter is None or not adapter.is_file():report['reason']='Archive identified, but verified default-EHT adapter unavailable';code=2
   else:
    report['adapter_sha256']=sha(adapter)
    for case in cases:
     # Send inputs only, never analytic eigenvalues or oracle labels/order.
     payload={k:case[k] for k in ['case','geometry_angstrom','basis','edges','electrons']};payload['archive_path']=str(archive.resolve());payload['route']='extended_huckel';payload['pyscf_allowed']=False
     # [K9] Adapter failures report exit code and stderr instead of a bare CalledProcessError.
     process=subprocess.run([sys.executable,str(adapter.resolve())],input=json.dumps(payload),text=True,capture_output=True,timeout=180)
     require(process.returncode==0,f"adapter failed on {case['case']} (exit {process.returncode}): {process.stderr.strip()[-2000:]}")
     result=json.loads(process.stdout);report['comparisons'].append(compare(case,result))
    missing=set(report['required_cases'])-{r['case'] for r in report['comparisons']}
    report['missing_cases']=sorted(missing)
    disagreements=[r['case'] for r in report['comparisons'] if r['status']!='PASS']
    report['status']='DISAGREEMENT_ORACLES_PREVAIL' if disagreements else ('BLOCKED' if missing else 'PASS')
    report['reason']='Missing release diatomic fixtures' if missing else ('Qualitative disagreements logged' if disagreements else 'All required comparisons computed')
    code=1 if disagreements else (2 if missing else 0)
  except Exception as e:report['status']='BLOCKED';report['reason']=f'{type(e).__name__}: {e}';code=2
 write_json(out/'checks/crystod_manifest.json',report);return code
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--archive',type=Path);p.add_argument('--adapter',type=Path);p.add_argument('--cases',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();code=run(a.out,a.archive,a.adapter,a.cases);print(json.loads((a.out/'checks/crystod_manifest.json').read_text())['reason']);sys.exit(code)
