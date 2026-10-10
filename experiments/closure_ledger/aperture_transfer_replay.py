"""Pinned replay plus disclosed post-measurement receiver-label correction.

The frozen producer's raw verdict remains NUMERICALLY_UNRESOLVED. This
module never alters its sources, archives, thresholds, or saved verdict.
"""
import copy
import json
from geometrodynamics.transaction import aperture_transfer as a
from experiments.closure_ledger import aperture_transfer_probe as p
import numpy as np

MANIFEST_SHA='26622f2692e5b53daa26c7207cde8548fea08874e4895d431d00b0787fd454d4'
SOURCE_ATOL=1e-14  # packet's analytic amplitude bound is one


def portable_record(record):
    """Validate the stored source, then normalize only transcendental roundoff.

    Authentication happens before this adapter in replay(). The original
    arrays are never mutated. Grid, inactive port and compact support remain
    exact checks; no relative tolerance is applied near packet zeros.
    """
    c=a.Config(**record['config']);a.validate(c)
    n=round((c.stop+.25)/c.dt)
    t=-.25+(np.arange(n)+.5)*c.dt
    if not np.array_equal(t,record['time']):raise ValueError('time grid mismatch')
    stored=np.asarray(record['incoming'])
    if stored.shape!=(n,2) or not np.all(np.isfinite(stored)):
        raise ValueError('invalid source array')
    local=np.zeros((n,2));local[:,c.source]=a.packet(t,c.carrier)
    if np.any(stored[:,1-c.source]!=0) or np.any(stored[abs(t)>=.25,c.source]!=0):
        raise ValueError('source support mismatch')
    error=float(np.max(abs(stored-local)))
    if error>SOURCE_ATOL:raise ValueError('source mismatch beyond roundoff tolerance')
    adapted=dict(record,incoming=local)
    return adapted,error


def portable_diagnose(record):
    return a.diagnose(portable_record(record)[0])


def authenticated_records(directory):
    """Keep the frozen byte-level authentication before numerical adaptation."""
    if p.digest(directory/'manifest.json')!=MANIFEST_SHA:raise ValueError('manifest mismatch')
    manifest=json.loads((directory/'manifest.json').read_text())
    names=[n+'.npz.b64' for n,c in p.schedule()]+['provenance.json','result.json']
    if set(manifest)!=set(names):raise ValueError('inventory mismatch')
    for n in names:
        if p.digest(directory/n)!=manifest[n]:raise ValueError('archive mismatch')
    provenance=json.loads((directory/'provenance.json').read_text())
    if provenance['sources']!=p.sources():raise ValueError('source mismatch')
    return {n:p.read(directory/(n+'.npz.b64')) for n,c in p.schedule()}


def reverse_receiver_diagnostic(record):
    """Use exact port-exchange covariance to diagnose arrival at A.

    In the real Gram factor, each row is (s,+s) or (s,-s); swapping ports
    therefore changes antisymmetric modal coordinates by -1. This checks
    the archived reversed-source trajectory, not a replacement trajectory.
    """
    c=a.Config(**record['config'])
    if c.source!=1 or not c.connected:raise ValueError('expected connected reverse-source control')
    _,W,_=a.operators(c)
    signs=np.where(W[:,0]==W[:,1],1.,-1.)
    r=copy.deepcopy(record);r['config']['source']=0
    r['incoming']=r['incoming'][:,::-1].copy();r['outgoing']=r['outgoing'][:,::-1].copy()
    r['final_q']*=signs;r['final_v']*=signs
    d=portable_diagnose(r)
    # The reversed experiment tests reciprocal bulk transport, not a reversed
    # MTY preparation. Do not report the relabelled virtual handle quantities.
    d.pop('advanced_return_fraction');d.pop('causal_return_fraction')
    d['physical_receiver']='A'
    return d


def replay(directory=p.RUN):
    records=authenticated_records(directory)
    normalized={};source_errors={}
    for name,record in records.items():
        normalized[name],source_errors[name]=portable_record(record)
    original=p.assess(normalized)
    saved=json.loads((directory/'result.json').read_text())
    if original['label']!=saved['label'] or original['checks']!=saved['checks']:
        raise ValueError('decision mismatch')
    corrected=copy.deepcopy(original)
    corrected['diagnostics']['reverse_source']=reverse_receiver_diagnostic(records['reverse_source'])
    # Reassess precisely the original numerical gates after replacing only
    # the source-relative arrival diagnostic. No tolerance changes.
    numeric=all(d['valid'] for d in corrected['diagnostics'].values())
    numeric=numeric and all(v<(1e-7 if n=='reverse_source' else .03) for n,v in corrected['comparisons'].items())
    numeric=numeric and corrected['unitarity_error']<1e-10 and corrected['reciprocity_error']<1e-10
    corrected['checks']['numerical_validity']=bool(numeric)
    corrected['label']=('FINITE_APERTURE_TRANSFER_SUPPORTED_AFTER_DIAGNOSTIC_CORRECTION'
                        if all(corrected['checks'].values()) else
                        'FINITE_APERTURE_TRANSFER_FAILED_AFTER_DIAGNOSTIC_CORRECTION'
                        if numeric else 'NUMERICALLY_UNRESOLVED')
    return dict(original_label=original['label'],original_checks=original['checks'],
                source_portability=dict(absolute_tolerance=SOURCE_ATOL,
                                        maximum_difference=max(source_errors.values()),
                                        differences=source_errors),
                correction='Post-measurement: reverse-source receiving port is A, not driven port B. No new trajectories or thresholds.',
                audited=corrected)


if __name__=='__main__':
    r=replay();print(json.dumps(dict(original_label=r['original_label'],correction=r['correction'],
                                   source_portability=r['source_portability'],
                                   audited_label=r['audited']['label'],checks=r['audited']['checks']),indent=2))
