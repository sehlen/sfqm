# Verification of the paper's computational tables

Target: the bundled August 18, 2014 extended *Lattices with many Borcherds
products*. This work builds on PR #4; master is not modified.

## Reproduce

Install the package with Sage 10.10, then run from outside the checkout:

```sh
sage -pip install --no-build-isolation /path/to/sfqm
cd /tmp
sage -python /path/to/sfqm/verification/verify_tables.py --entries 9 --output table9.json
sage -python /path/to/sfqm/verification/verify_tables.py --search 3 --output closure_n3.json
```

Entry jobs cover Tables 4 (75), 5 (79 weight/symbol pairs), 7 (70), 8 (43),
9 (169), and the 90 cyclic modules listed on page 25. The manifests include
signatures, rank bounds, and expected global realizability where applicable.
The exact dimension expression uses QQbar and refuses to round a noninteger.
Weight 2 and 3/2 corrections first use good-reduction upper bounds to certify
zero; inconclusive cases fall back to characteristic-zero invariant computations.
The character/value routines and invariant implementation are shared with the
library, so this is not an independent implementation of all underlying math.

The extension searches implement Algorithm 5.5 from the published anisotropic
seeds, with every prime below Table 6's bound and without dynamic prime pruning.
They compare actual isomorphism classes, not just counts. A completed closure
is conditional on the published seed list being complete. The optional
`--seed-level` scan checks an additional bounded range, not the full theoretical
bounds. `complete` in the JSON means that job finished, not an unconditional
proof of the classification. Full exclusion of higher-level anisotropic seeds,
Table 1's genus-to-isometry multiplicities, and the structural Tables 2, 3, 6
require separate verification. Do not claim the 362-lattice theorem certified
from these jobs alone.

JSON checkpoints are written after each row / completed graph vertex, and
uploaded as GitHub Actions artifacts even on failure. Logs include exact rows
and a final SUMMARY. Timeouts or runtime errors are failures, not passing checks.

## Data provenance

`extract_manifests.py /path/to/dimest.tex` recreates the manifests using the
arXiv v1 TeX source (https://arxiv.org/src/1408.4148v1) and the bundled extended
PDF. It requires PyMuPDF, not Sage. Table 8 was transcribed explicitly; Table 9
uses PDF glyph positions to distinguish superscripts and subscripts. Input hashes
are in `sources.json`. Tables 4, 5 and 7 were extracted from the v1 TeX and
correspond to the tables read in the bundled PDF. Extraction assertions check
row counts; mathematical verification separately checks each entry.

The first local attempt checked Table 8 and reproduced its n=3 closure, but
workspace maintenance removed those unpublished result files. The durable CI
runs are the recorded evidence for this verification branch.

## Why reduction can certify zero

Let b be the exactly computed virtual cusp dimension, and c the nonnegative
low-weight correction. The desired dimension is b+c. After reduction at a
splitting prime coprime to the level and group order, the rank of the invariant
equations cannot increase, so the computed correction u satisfies c <= u.
The +/- eigenspaces in the weight-one-half correction obey the same bound
since the residue characteristic is odd. Consequently b+u=0, together with
nonnegativity of the true dimension, proves b+c=0. This is a zero certificate,
not a general assertion that modular and characteristic-zero kernels agree.
When b+u>0 the verifier computes the exact characteristic-zero correction.
A negative b+u raises an error. Each entry records which method certified it.

The original all-characteristic-zero run remains available separately; the
follow-up workflow targets the slower weight-3/2 jobs. Full entry and closure
commands for other tables remain available through the script.
