import json
import numpy as np
from jsonschema import Draft202012Validator
from common import csv_read
class Parameters:
 def __init__(self,root):
  self.root=root;self.on={}
  for r in csv_read(root/'params/onsite.csv'):self.on[(r['element'],r['orbital'][-1])]=float(r['energy_eV'])
  self.rows=csv_read(root/'params/pairs.csv');self.defaults=json.loads((root/'params/channel_defaults.json').read_text())['defaults']
  data=json.loads((root/'params/presets.json').read_text());schema=json.loads((root/'params/presets.schema.json').read_text());Draft202012Validator(schema).validate(data)
  self.presets={p['id']:p for p in data['presets']};assert len(self.presets)==len(data['presets'])
  seen=set()
  for r in self.rows:
   assert r['distance_law'] in ['harrison_inverse_square','fixed','off'];assert r['provenance_kind'] in ['harrison','harrison_extrapolated_H','pedagogical_override','lecture_decoupling'];assert r['s_carrier'] in ['H','none']
   key=tuple(sorted([r['element_a'],r['element_b']]))+(r['channel'],);assert key not in seen;seen.add(key)
   assert r['source'] and r['units']=='eV'
   if r['provenance_kind'] in ['pedagogical_override','lecture_decoupling']:assert r['override_reason']
   if r['distance_law']=='fixed':float(r['value_eV'])
   if r['distance_law']=='harrison_inverse_square':float(r['eta']);assert float(r['K_eV_A2'])>0
   assert float(r['scale'])>0
 def onsite(self,e,o):return self.on[(e,o)]
 def get(self,a,b,ch,d,preset):
  p=self.presets[preset];rows=[r for r in self.rows if sorted([r['element_a'],r['element_b']])==sorted([a,b])];hits=[r for r in rows if r['channel']==ch]
  # [K22] CONVENTIONS v1.6 W2/W3: universal-law default only for a pair with no row at all, both elements with on-sites;
  # the H default when the pair contains H. A pair with rows keeps "missing channel is an error".
  if not rows and a in {e for e,_ in self.on} and b in {e for e,_ in self.on}:hits=[r for r in self.defaults if r['channel']==ch and r['s_carrier']==('H' if 'H' in [a,b] else 'none')]
  if len(hits)!=1:raise ValueError(f'No unique pair/default row for {a}-{b} {ch}; absent is not zero')
  r=hits[0];law=r['distance_law'];scale=float(r['scale'])
  if r['s_carrier']=='H':
   assert 'H' in [a,b];assert abs(scale-p['hydrogen_scale'])<1e-15,'preset/row H scale conflict'
  if law=='off':return 0.
  if law=='fixed':return float(r['value_eV'])*scale
  if d<=0:raise ValueError('nonpositive distance')
  return float(r['eta'])*float(r['K_eV_A2'])*scale/d**2
 def element(self,a,oa,b,ob,u,d,preset):
  for ov in self.presets[preset]['coupling_overrides']:
   direct=(a,oa,b,ob)==(ov['element_a'],ov['orbital_a'],ov['element_b'],ov['orbital_b']);rev=(b,ob,a,oa)==(ov['element_a'],ov['orbital_a'],ov['element_b'],ov['orbital_b'])
   if direct or rev:return float(ov['value_eV'])
  if oa=='s' and ob=='s':return self.get(a,b,'ss_sigma',d,preset)
  if oa=='s':return self.get(a,b,'sp_sigma',d,preset)*u['xyz'.index(ob)]
  if ob=='s':return -self.get(a,b,'sp_sigma',d,preset)*u['xyz'.index(oa)]
  i='xyz'.index(oa);j='xyz'.index(ob);ps=self.get(a,b,'pp_sigma',d,preset);pp=self.get(a,b,'pp_pi',d,preset)
  return (pp if i==j else 0.)+(ps-pp)*u[i]*u[j]
 def block(self,a,b,u,d,preset):
  aa=['s'] if a=='H' else ['s','x','y','z'];bb=['s'] if b=='H' else ['s','x','y','z']
  return np.array([[self.element(a,x,b,y,u,d,preset) for y in bb] for x in aa])
