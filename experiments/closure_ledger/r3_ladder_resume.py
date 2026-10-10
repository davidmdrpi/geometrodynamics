"""Resume the ladder run after the 2-hour background-command limit stopped it (docs/r3_ladder.md).

The first production process (started 18:27 UTC) was killed at 20:27 UTC during the 4/11 scan; the six
scans already written are kept as written. This harness calls only the frozen r3_ladder_probe functions
(_scan_one, _hp_one, _hp2_one, bracket, _write, utc) and reproduces each frozen per-rung loop body and
record format exactly, one rung at a time, skipping any archive that already exists, so that each call
fits within the limit. It changes no setting, tolerance, rung or phase.
Usage: python -m experiments.closure_ledger.r3_ladder_resume LABEL STAGE [STAGE ...]   (STAGE: scan hp hpnoise score)
"""
import sys
from multiprocessing import Pool
from experiments.closure_ledger import r3_ladder_probe as lp


def run(stage):
    if stage == 'score':
        if not (lp.RUN_DIR/'result.json').exists():
            lp.main(['score'])
        return
    name, fn = dict(scan=('scan', lp._scan_one), hp=('hp', lp._hp_one), hpnoise=('hpnoise', lp._hp2_one))[stage]
    for p, q in lp.RUNGS:
        t = lp.tag(p, q)
        if (lp.RUN_DIR/f'{name}_{t}.json').exists():
            continue
        start = lp.utc()
        with Pool(4) as pool:
            out = pool.map(fn, [(p, q, j) for j in range(lp.N_PHASE)])
        if stage == 'scan':
            br = lp.bracket(p, q)
            rec = dict(p=p, q=q, started_utc=start, finished_utc=lp.utc(), points=out,
                       **{k: v for k, v in br.items() if k != 'K'})
        else:
            rec = dict(p=p, q=q, started_utc=start, finished_utc=lp.utc(), rows=out)
        lp._write(f'{name}_{t}.json', rec)
        print(f'{name} {p}/{q} done', flush=True)


def main(argv):
    label, stages = argv[0], argv[1:]
    lp._write(f'resume_{label}.json', dict(started_utc=lp.utc(), stages=stages,
                                           note='resumption after background time limit; frozen functions only'))
    for s in stages:
        run(s)
        print(f'stage {s} finished', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
