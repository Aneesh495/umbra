# ADR 0008: Cryptographic Evidence Bundling and Adversarial Verification

## Status
Accepted

## Context
Scientific computational codes must guarantee reproducibility, numerical fidelity, and provenance. Without automated verification mechanisms, silent numerical regressions, uncalibrated parameter shifts, or unauthorized source modifications can compromise simulation integrity unnoticed.

To establish verifiable confidence in `umbra`, we required an automated verification protocol that links numerical benchmark output directly to cryptographic signatures of the exact source tree.

## Decision
We implement a cryptographically signed verification architecture comprising:

1. **Automated Evidence Generation (`scripts/run_acceptance.py`)**:
   Runs the full suite of acceptance benchmarks:
   - Relativistic Sod Shock Tube
   - Relativistic Brio-Wu Magnetized Shock
   - Relativistic Orszag-Tang Vortex
   - Kerr Black Hole Geodesic Raytracer
   - Fishbone-Moncrief Equilibrium Torus
2. **Deterministic Metric Collection**:
   Extracts strict quantitative convergence metrics:
   - Maximum magnetic divergence $\max |\nabla \cdot \mathbf{B}|$
   - Conservation error of Carter constant $\Delta Q / Q_0$
   - Black hole shadow radius $R_{shadow}$ compared against the analytical Bardeen formula
   - L1 and L2 norm convergence rates across shock and smooth flows
3. **Source Tree SHA-256 Fingerprinting**:
   Calculates a single SHA-256 digest encompassing all Fortran, C, and Python source files.
4. **Cryptographic Payload Signing (`artifacts/evidence_bundle.json`)**:
   Serializes all metrics, environment fingerprints (compiler version, architecture, timestamp), and source digests into a canonical JSON payload, hashed with SHA-256.
5. **Adversarial Negative Control Suite (`scripts/verify_evidence.py` & `tests/negative/`)**:
   Executes 12 independent negative controls to audit:
   - Rejection of tampered evidence payloads
   - Detection of illegal syntax or style infractions (such as em dash presence)
   - Enforcement of unphysical super-extremal spin rejection ($a > 1$)
   - Preservation of physical density and pressure floors
   - Bit-for-bit matching of live repository files against bundle fingerprints

## Consequences
- **Advantages**:
  - Independent third parties can audit simulation results and verify source tree authenticity with a single command (`python3 scripts/verify_evidence.py`).
  - Guarantees zero regressions across compiler versions or hardware platforms.
  - Hardens the codebase against subtle mathematical and physical boundary condition errors.
- **Trade-offs**:
  - Requires updating the evidence bundle when intentional algorithmic improvements alter the source tree.
