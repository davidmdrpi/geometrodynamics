"""Pinned diagnostic replay of the immutable finite-packet MTY experiment."""
import json
from experiments.closure_ledger import mty_packet_probe as p

MANIFEST_SHA='c2156a1de4ed711204f90ac333e8dfc933b0dd4223cb2e1008613c9810fcca2d'


def replay(directory=p.RUN):
    return p.replay(directory,MANIFEST_SHA)


if __name__=='__main__':
    r=replay()
    print(json.dumps({k:r[k] for k in ('label','checks','gr_traversable_support',
          'gravitating_mouth_recoil','unique_history','action_discreteness')},indent=2))
