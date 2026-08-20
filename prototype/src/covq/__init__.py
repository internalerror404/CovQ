"""CovQ: quantum Fisher-information compilation for commuting Pauli programs.

Pre-pilot research code.  Three conventions hold everywhere and are repeated in
the modules that depend on them:

* QFIMs are **per shot**;
* resource counts come from **emitted, lowered (and routed) circuits**, never
  from a formula;
* the correctness gates are separated from the scientific hypotheses, and no
  tolerance is adjusted in response to a failure.
"""

from . import circuits, gates, instances, paulis, polytope, programs, qfim, sim, topology, width

__all__ = [
    "circuits", "gates", "instances", "paulis", "polytope", "programs",
    "qfim", "sim", "topology", "width",
]
__version__ = "0.1.0"
