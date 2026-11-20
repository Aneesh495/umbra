# ADR 0008: Source-Bound Evidence and Verification

## Status

Accepted

## Context

The verification runner needs to record which source produced a result and the measurements produced by that run. A saved evidence bundle also needs a way to detect accidental changes to its contents.

## Decision

`scripts/run_acceptance.py` executes the unit and acceptance binaries and records their measured output. The five acceptance checks cover Sod and Brio-Wu initial states, initial magnetic divergence in the Orszag-Tang vortex, Carter constant drift, and Schwarzschild shadow radius. The bundle records a SHA-256 digest of the supported files under `src/` and a SHA-256 digest of the serialized payload.

`scripts/verify_evidence.py` checks those digests, the tests' reported counts and tolerances, and the local negative controls. The `cryptographic_signature` field is the historical schema name for the payload digest. Publisher authentication would require a separate signature and a trusted verification key.

## Consequences

`make test` builds the required executables, generates a temporary bundle from the current sources, and runs the verifier. The checked-in example bundle records its original revision and measurements. To retain a fresh bundle separately, run `python3 scripts/run_acceptance.py --output path/to/evidence.json` after building the test binaries.
