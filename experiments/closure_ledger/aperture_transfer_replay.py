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
    d=a.diagnose(r)
    # The reversed experiment tests reciprocal bulk transport, not a reversed
    # MTY preparation. Do not report the relabelled virtual handle quantities.
    d.pop('advanced_return_fraction');d.pop('causal_return_fraction')
    d['physical_receiver']='A'
    return d


def replay(directory=p.RUN):
    original=p.replay(directory,MANIFEST_SHA)
    corrected=copy.deepcopy(original)
    corrected['diagnostics']['reverse_source']=reverse_receiver_diagnostic(p.read(directory/'reverse_source.npz.b64'))
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
                correction='Post-measurement: reverse-source receiving port is A, not driven port B. No new trajectories or thresholds.',
                audited=corrected)


if __name__=='__main__':
    r=replay();print(json.dumps(dict(original_label=r['original_label'],correction=r['correction'],
                                   audited_label=r['audited']['label'],checks=r['audited']['checks']),indent=2))
