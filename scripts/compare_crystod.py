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

def compare(case,result):
 assert result['route']=='extended_huckel' and result.get('pyscf_used') is False
 assert result['case']==case['case']
 Q=np.array(result['axis_map_checker_to_contract'],float);assert Q.shape==(3,3) and np.max(abs(Q.T@Q-np.eye(3)))<1e-10
 assert result['basis_description'] and result['basis_differences']
 # Explicit label map must be backed by the adapter's actual coordinate convention.
 mapping=result['irrep_map_to_contract'];assert isinstance(mapping,dict) and mapping
 got=result['levels'];assert all(x['irrep'] in mapping for x in got)
 expected=case['levels'];shared=set(mapping[x['irrep']] for x in got)&{x['irrep'] for x in expected};assert shared
 full_expected='Dinfh' if case['case']=='H2' else case['group'];point_ok=result['point_group_full']==full_expected
 checks=[]
 for label in sorted(shared):
  ours=[x for x in expected if x['irrep']==label];theirs=sorted([x for x in got if mapping[x['irrep']]==label],key=lambda x:x['energy'])
  # canonical_level_id is independently assigned by the adapter (branch/correlation identity), not inferred from the ordering being tested.
  ours_ids=[f'{label}:{k+1}' for k in range(len(ours))]
  ids=[x['canonical_level_id'] for x in theirs]
  dims=[x['multiplicity'] for x in theirs];want=[x['multiplicity'] for x in ours]
  checks.append({'irrep':label,'shared_basis_counts_match':len(ours)==len(theirs),'degeneracy_pattern_matches':dims==want,'label_after_axis_mapping':label,'order_within_irrep_matches':ids==ours_ids,'expected_order':ours_ids,'checker_order':ids,'identity_assignment_evidence':result['level_identity_evidence']})
 ok=point_ok and all(c['shared_basis_counts_match'] and c['degeneracy_pattern_matches'] and c['order_within_irrep_matches'] for c in checks)
 return {'case':case['case'],'geometry_angstrom':case['geometry_angstrom'],'contract_basis':case['basis'],'checker_basis':result['basis_description'],'basis_differences':result['basis_differences'],'axis_map_checker_to_contract':Q.tolist(),'irrep_map':mapping,'point_group_matches':point_ok,'expected_full_point_group':full_expected,'checker_full_point_group':result['point_group_full'],'shared_irrep_checks':checks,'status':'PASS' if ok else 'DISAGREEMENT_ORACLES_PREVAIL','energy_comparison':'NOT_PERFORMED; raw energies used only to sort within each irrep','composition_and_bond_orders':'NOT_COMPARED'}

def run(out,archive=None,adapter=None,cases_file=None):
 cases=json.loads((cases_file or ROOT/'checks/crystod_cases.json').read_text())
 report={'schema_version':'1.2','status':'BLOCKED','archive_filename':None,'archive_sha256':None,'version':None,'comparisons':[],'allowed':['point_group','shared_irrep_degeneracy','axis_mapped_labels','within_irrep_order'],'prohibited':['energy_values','compositions','bond_orders'],'precedence':'Appendix A oracles win; discrepancies recorded for the owner','engine_UI_gate':'separate, not exercised by this script','diatomic_gate':'B2,C2,O2,F2,CO missing; N2 scope undecided','required_cases':[x['case'] for x in cases]+['B2','C2','O2','F2','CO']}
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
     process=subprocess.run([sys.executable,str(adapter.resolve())],input=json.dumps(payload),text=True,capture_output=True,timeout=180,check=True)
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
