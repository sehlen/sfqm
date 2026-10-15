# Paper-table verification results

The target is the repository's August 18, 2014 extended paper. Results are
recorded as machine-readable JSON with the Sage version, code commit, job/run
IDs, elapsed computation time, and entry-by-entry checks or missing/extra lists.
Master remains unchanged. This verification is proposed in PR #6 atop PR #4.

## Completed characteristic-zero checks

The first run is [38685267716](https://github.com/sehlen/sfqm/actions/runs/38685267716).

- Table 5: all 79 weight/module pairs pass, including signature and anisotropy.
- Table 7: all 70 entries pass, including the three exclusions from global realization.
- Table 8: all 43 entries pass, including the four exclusions from global realization.
- The n=2 extension search reproduces all 70 Table 7 classes, with 67 globally
  realizable genera. It takes about 342 seconds of computation on its runner.
- The completed higher-weight closures agree by isomorphism class, not merely
  by count; see the individual `results/closure_n*.json` files.

The weight-3/2 entry checks and closure are being rerun with good-reduction
zero certificates in [38686028123](https://github.com/sehlen/sfqm/actions/runs/38686028123).
The original characteristic-zero-only jobs are retained for comparison.
No uncompleted job is counted as a successful verification.

## Meaning and limits

Entry checks establish the printed modules' vanishing dimensions, signatures,
rank bounds and, where applicable, global existence/exclusion tests. They share
the library's character and value routines, so are not a second implementation
of every mathematical ingredient. Their dimension arithmetic avoids numerical
rounding. The low-weight certificate and its fallback are explained in README.md.

Closure implements the extension stage of Algorithm 5.5 from the published
anisotropic seed lists. Its completeness is conditional on those lists being
complete. The full theoretical-range exclusion of further anisotropic seeds has
not been rerun here. Structural Tables 2, 3 and 6 are used as theoretical input,
not independently proved by these computations.

Table 1 reconciles with the extracted genera: 67 for n=2; 15,5,3,4,2,3,3,2 for
n=3,...,10; one each for n=18 and 26. For n=1, the 165 realizable entries of
Table 9 plus 90 cyclic genera give 255 genera. Proposition 6.2 contributes one
extra class in the exceptional genus 2_1^+1.5^-1.25^-1, giving 256 and total 362.
This is an accounting check conditional on that proposition, not a new class
number calculation. `summarize.py` reproduces this accounting and the run summary.

For n>=2, the largest module orders in Tables 7 and 8 are 625 and 512,
respectively, consistent with the determinant threshold used in Proposition 6.2.

The divisor-sum sigma bug from PR #4 is distinct from the Gauss-sum character
used in the dimension formula. In the reviewed old code the sigma-based bound
was commented out; this defect alone did not alter that active dimension path.
