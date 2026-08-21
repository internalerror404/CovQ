"""Is the registered tE ratio on the deep instances tight?

`bE` could not be afforded on `path_m4`, `path_m5` or `banded_m5` -- exactly the
three largest emitted-CX ratios -- so their upper end rests on a single arm.
This probes that arm with cheaper alternatives.  It is a *falsification* attempt,
not a search for a better number to report: if anything beats the registered
run, the registered ratio is an overestimate and the manuscript must say "up to".

Arms
----
``tE`` from ``|+>^m``      the registered configuration.
``bE`` windowed at ``W``   best-position search restricted to the last ``W``
                           slots.  ``W = 1`` is exactly tE and ``W = inf`` is
                           exactly bE, so the *cost* interpolates -- linear in
                           depth rather than quadratic.  The *quality* does not
                           interpolate, which is the point of running it.
``tE`` random restarts     Haar-ish product-free starts, to test whether the
                           ``|+>^m`` rule is doing real work.

Findings this reproduces
------------------------
- ``W = 4`` beats the registered tE on ``path_m4``: 20 emitted CX against 24.
  The registered deep-instance ratios are therefore overestimates.
- ``W = 2`` is *worse* than ``W = 1`` on both instances where it ran, and on
  ``path_m4`` it does not even converge.  Greedy insertion is not monotone in
  its candidate set, so no windowed run may ever be quoted as a bound on bE.
- Random restarts are roughly four times worse than ``|+>^m`` (102 and 110 CX
  against 24), which independently justifies the registered initialisation.
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from covq.circuits import lower_to_cx, resources
from covq.instances import banded, path_target
from covq.paulis import z_generators
from covq.programs import labelled_schedule_from_width
from covq.quest import moment_constraints, quest_circuit, quest_published
from covq.width import decompose_width2

DEEP = {
    "path_m4": (lambda: path_target(4, 0.45).F, 4),
    "path_m5": (lambda: path_target(5, 0.45).F, 5),
    "banded_m5": (lambda: banded(5, 2, 0.3).F, 5),
}


def arms(m, windows, restarts, rng):
    out = [("tE (registered)", dict(variant="tE"))]
    out += [(f"bE window={w}", dict(variant="bE", position_window=w)) for w in windows]
    for k in range(restarts):
        v = rng.normal(size=1 << m) + 1j * rng.normal(size=1 << m)
        out.append((f"tE restart{k}", dict(variant="tE", psi0=v / np.linalg.norm(v))))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--instances", nargs="*", default=list(DEEP))
    ap.add_argument("--windows", nargs="*", type=int, default=[4])
    ap.add_argument("--restarts", type=int, default=0)
    ap.add_argument("--max-depth", type=int, default=80)
    args = ap.parse_args()

    print(f"{'instance':<12}{'arm':<20}{'rot':>5}{'cx':>5}{'ratio':>8}{'secs':>7}  verdict")
    for name in args.instances:
        make, m = DEEP[name]
        F = make()
        dec = decompose_width2(F)
        covq_cx = sum(w * sum(1 for b in blocks if len(b) == 2)
                      for w, blocks, _ in dec.branches)
        cons = moment_constraints(z_generators(m), F)
        best = None
        for tag, kw in arms(m, args.windows, args.restarts, np.random.default_rng(11)):
            t0 = time.time()
            r = quest_published(cons, m, max_depth=args.max_depth, tol=1e-13, **kw)
            cx = resources(lower_to_cx(quest_circuit(r, m)))["cx_count"]
            if tag.startswith("tE (reg"):
                best = cx
                verdict = "registered"
            elif r.stop_reason != "converged":
                verdict = f"no verdict ({r.stop_reason})"
            elif cx < best:
                verdict = f"BEATS registered by {best - cx} CX"
            else:
                verdict = "does not beat registered"
            print(f"{name:<12}{tag:<20}{r.depth_adaptive_length:>5}{cx:>5}"
                  f"{cx / covq_cx:>8.2f}{time.time() - t0:>7.0f}  {verdict}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
