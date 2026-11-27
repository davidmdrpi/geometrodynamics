"""Case 2 of the Bloch-sphere joint condition, in the weak-coupling limit (review argument, checked here).

Both visits free (t = 1/4, 3/4); the singlet lock sets e_B = -e_A = -e with e perpendicular to n0, chosen
as the leading eigenvector of a_perp a_perp^T + b_perp b_perp^T. Outcomes A = sgn(e.a), B = sgn(-e.b).
Reports CHSH, the binary marginals, and Alice's E[(2u_A - 1)^2] = E[(e.a)^2] as Bob's setting changes.
"""
import json
from pathlib import Path
import numpy as np

rng = np.random.default_rng(2027)
n = rng.normal(size=(400_000, 3)); n /= np.linalg.norm(n, axis=1, keepdims=True)


def axis(phi):
    return np.array([np.cos(phi), np.sin(phi), 0.0])


def lock(a, b):
    ap = a - (n @ a)[:, None] * n
    bp = b - (n @ b)[:, None] * n
    M = ap[:, :, None] * ap[:, None, :] + bp[:, :, None] * bp[:, None, :]
    w, v = np.linalg.eigh(M)
    e = v[:, :, -1]                                  # leading eigenvector, orthogonal to n0 (n0 is in the kernel)
    e *= np.where(rng.random(len(n)) < 0.5, 1, -1)[:, None]   # K_A x K_B pairs e with -e
    return e


def stats(a, b):
    e = lock(a, b)
    A, B = np.sign(e @ a), np.sign(-(e @ b))
    return dict(E=float(np.mean(A * B)), pA=float(np.mean(A > 0)), pB=float(np.mean(B > 0)),
                alice_intensity=float(np.mean((e @ a) ** 2)))


a0, a1, b0, b1 = axis(0.0), axis(np.pi / 2), axis(np.pi / 4), axis(-np.pi / 4)   # Bloch-sphere CHSH axes
pairs = {"a0b0": stats(a0, b0), "a0b1": stats(a0, b1), "a1b0": stats(a1, b0), "a1b1": stats(a1, b1)}
S = pairs["a0b0"]["E"] + pairs["a0b1"]["E"] + pairs["a1b0"]["E"] - pairs["a1b1"]["E"]
signal = dict(parallel=stats(a0, a0)["alice_intensity"], perpendicular=stats(a0, a1)["alice_intensity"])
out = dict(samples=len(n), pairs=pairs, S=float(S), alice_intensity_vs_bob=signal)
Path(__file__).with_name("bloch_joint_case2.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
