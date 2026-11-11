# How many qubits and gates does Shor's algorithm need?

A small, rerunnable resource estimate for Shor's algorithm. It goes in three steps:

1. **Known answer.** Factor N = 15 (and N = 21) in Qiskit with every modular multiplication built
   from real adders, count the compiled gates, simulate the circuit and recover the factors. Every
   number in this step is measured.
2. **Scale up.** Estimate logical qubits, Toffoli gates, code distance, physical qubits and runtime
   for an n-bit number, up to RSA-2048, with a simple model whose assumptions are written down below.
3. **Check the rules.** Build and simulate small cases to test each rule the scale-up relies on.

```bash
pip install -r requirements.txt                      # Python 3.9+, versions pinned
python3 shor_arithmetic.py --json arith.json         # N = 15 from real arithmetic (about 40 s)
python3 shor_arithmetic.py --N 21 --a 2 --json arith21.json   # N = 21, 20 qubits (about 2-3 min)
python3 shor_estimate.py --json results.json         # comparison oracle for N = 15 + table for 8 ... 4096 bits
python3 shor_estimate.py --skip-a --bits 2048 --p 2e-3        # change the size or the gate error
python3 verify_rules.py --json rules.json            # rule checks (about 35 s)
python3 figures/make_shor_rules.py out               # any figure script; writes PNG + JPG to out/
```

Run the scripts from this folder: `shor_arithmetic.py` and `verify_rules.py` import `shor_estimate.py`.
The `*.json` files are the saved outputs of these commands; the figure scripts in `figures/` read them.

The five social-media figures are in [`figures/images/`](figures/images/); the scripts next to them rebuild them from the saved `*.json` outputs.

## 1. What is measured, what is taken from papers, what is estimated

| Origin | Numbers | Trust |
|---|---|---|
| **Measured** (this repo, Qiskit) | N = 15 from real arithmetic: 18 qubits, 14,256 CX, depth 21,689, peaks at 0/64/128/192, r = 4, 15 = 3 × 5. One controlled ×7 mod 15: 1,784 CX. N = 21: 20 qubits, r = 6, 21 = 7 × 3 | Exact for this circuit |
| **From papers** (formulas we do not re-derive) | logical qubits 3n + 0.002 n lg n, Toffolis 0.3n³ + 0.0005 n³ lg n, measurement depth 500n² + n² lg n, logical-error rule of thumb, 10 µs reaction time [GE19]; minimum logical qubits 2n + 3 [B03] | Peer-reviewed design |
| **Estimated** (our model on top) | code distance d, physical qubits, hours | Order of magnitude only |

`shor_arithmetic.py` builds every controlled multiplication from Draper QFT adders, a modular adder
and a modular multiplier (Beauregard 2003), for any a and N. It checks the controlled multiplier on
every input x < N with the control on and off, keeps all controlled multiplications (even when
a^(2^j) mod N = 1), and validates the period (r even, a^(r/2) ≠ −1 mod N).

`shor_estimate.py` Part A builds the same algorithm with a **precomputed-permutation oracle**: each
controlled multiplication is a permutation matrix built from the full classical multiplication
table, which is affordable for N = 15 and impossible for 2048-bit numbers. It is not given the
period, but the compiler drops the multiplications that equal ×1, and those follow from it. Like for like, one controlled ×7
mod 15 costs **1,784 CX** with real arithmetic and **236 CX** as a permutation (7.6×), both at
optimization level 1. Whole circuits: 18 qubits and 14,256 CX against 12 qubits and 540 CX; most
of that gap comes from the six ×1 multiplications the oracle version drops. Small demos that go
further and simplify the circuit with knowledge of the period are criticised in Smolin, Smith,
Vargo, "Pretending to factor large numbers on a quantum computer", arXiv:1301.7007.

## 2. The spec sheet of our model, and why

A resource estimate means nothing without the machine it targets. Machine specs come from [GE19];
three lines are **our model choices** and are marked as such.

| Spec | Value | Why this value | If it gets worse |
|---|---|---|---|
| Physical error rate p (every gate, measurement, idle step) | **10⁻³** (0.1 %), optimistic case 10⁻⁴ | The assumption in [GE19] and [G25], so results are comparable. | d grows and qubits grow fast (table below). |
| Errors | independent, uniform, no leakage, no correlated bursts | The textbook surface-code model. | Correlated errors can break the scaling law. |
| Code | rotated surface code, 2D nearest-neighbour grid, lattice surgery | High threshold with only local 2D links, which matches superconducting chips. | qLDPC codes can need fewer qubits but need longer-range links. |
| Logical error per qubit per round | **0.1 · (p / 1 %)^((d+1)/2)** | Rule of thumb of [GE19] (Sec. 2.13, Fig. 8). The 1 % is a property of code + noise + decoder, not of the machine. Each step d → d + 2 divides the error by Λ = 1 % / p, so **Λ = 10 at p = 10⁻³**, as [GE19] assumes. | At p = 0.5 % (Λ = 2 in this model) d = 101 and about 190 million qubits. |
| Qubits per logical qubit | **2d² − 1 ≈ 2d²** (d² data + d² − 1 measure qubits); [GE19] uses 2(d + 1)² | Rotated surface-code layout, checked in `verify_rules.py`. | |
| Routing space (model choice) | **× 1.5** on top of the data patches | Lattice surgery needs free patches to move and connect logical qubits. | Linear in qubit count. |
| Magic-state factories (model choice) | **4 factories × 12 patches** | Crude; [GE19] and [G25] model factories in detail. In the [GE19] layout they take about a quarter of the chip. | Underestimates factories (see limitations). |
| Round time (one round of stabilizer measurements) | **1 µs** | Superconducting rounds are about 1 µs; [GE19] uses the same. | Runtime scales linearly: 1 ms rounds (reaction scaled to 10 rounds = 10 ms) turn 6 h into about 8 months; with reaction kept at 10 µs, about 25 days. |
| Reaction time (decode a measurement and choose the next step) | **10 µs = 10 rounds per logical layer** | [GE19] is reaction-limited: about 10 µs per layer of its measurement depth. Needs a decoder that keeps up in real time. | A slow decoder stalls the computation. |
| Failure budget (model choice) | **1 %** chance that the whole run fails | Strict on purpose. [GE19] accepts 31 % per run and retries; [G25] has 6.7 % per shot and expects 9.2 shots. | Budget 0.1 % adds about one step of d. |

### Measured error suppression so far

Λ is the factor by which the logical error drops each time the code distance grows by 2. These
values come from different codes, protocols and decoders, so they are not a league table.

| Group, platform | Code, protocol | Λ | Source | Caveat |
|---|---|---|---|---|
| USTC, superconducting (Zuchongzhi 3.2) | surface code, d up to 7 | 1.40 ± 0.06 | He et al., PRL 135, 260601 (2025) | |
| Google, superconducting (Willow, 105 qubits) | surface code, d up to 7 | 2.14 ± 0.02 | Google Quantum AI, Nature 638, 920 (2025) | 1.1 µs rounds; real-time decoding shown at d = 5 |
| Harvard/MIT/QuEra, neutral atoms (448 atoms) | surface code, d = 3 vs d = 5, 4-round circuit | 2.14 ± 0.13 | Bluvstein et al., Nature (2025), arXiv:2506.20661 | a two-point ratio, not a multi-distance memory fit |
| Quantinuum, trapped ions (Helios) | **color code**, d up to 7 | 12 (+3/−2) | arXiv:2610.24077 (2026, preprint) | different code; 73 qubits; Λ from offline decoding (real-time-capable decoders benchmarked); round time not stated (one leakage-repump cycle alone takes about 340 µs) |

The estimates need surface-code Λ around 10 **and** 1 µs rounds with real-time decoding **and**
about a million qubits at the same time. No platform shows all three yet: the superconducting
surface-code chips have the speed but Λ ≈ 1.4–2.1; the ion color-code result reaches a similar
suppression on a different code, but not yet at that speed or scale.

## 3. How the estimate is computed

For an n-bit number N:

1. Logical qubits Q = 3n + 0.002 n lg n, Toffolis T = 0.3n³ + 0.0005 n³ lg n, and logical layers
   D = 500n² + n² lg n, all from [GE19].
2. Each layer waits one reaction time (10 rounds), so the run takes D · 10 rounds.
3. Pick the smallest odd d with Q · D · 10 · p_L(d) < 1 %: "number of qubit-rounds times error per qubit-round".
4. Physical qubits = Q · 2d² · 1.5 + 4 · 12 · 2d².
5. Hours = D · 10 · 1 µs.

Example, RSA-2048 at p = 10⁻³: Q = 6,189 and D = 2.14·10⁹, so p_L must stay below
7.5·10⁻¹⁷ per qubit per round. d = 29 gives 10⁻¹⁶ (too high) and d = 31 gives 10⁻¹⁷, so
**d = 31**. One logical qubit then costs 2 · 31² = 1,922 physical qubits, or about 2,900 with the
×1.5 routing space ([GE19]: 2 · 28² = 1,568 at d = 27, routing counted separately in its layout).

## 4. Results

| Bits | Logical qubits | Toffolis | d (p = 10⁻³) | Physical qubits | Hours |
|---|---|---|---|---|---|
| 256 | 772 | 5.1·10⁶ | 25 | 1.5·10⁶ | 0.09 |
| 1024 | 3,092 | 3.3·10⁸ | 29 | 7.9·10⁶ | 1.5 |
| 2048 | 6,189 | 2.6·10⁹ | 31 | 1.8·10⁷ | 6.0 |
| 4096 | 12,386 | 2.1·10¹⁰ | 33 | 4.1·10⁷ | 24 |

RSA-2048 versus the gate error (same model; runtime does not depend on d because layers are reaction-limited):

| p | Λ in this model | d | Physical qubits | Hours |
|---|---|---|---|---|
| 10⁻⁴ | 100 | 15 | 4.2·10⁶ | 6.0 |
| 10⁻³ | 10 | 31 | 1.8·10⁷ | 6.0 |
| 2·10⁻³ | 5 | 43 | 3.5·10⁷ | 6.0 |
| 3·10⁻³ | 3.3 | 57 | 6.1·10⁷ | 6.0 |
| 5·10⁻³ | 2 | 101 | 1.9·10⁸ | 6.0 |

**Comparison with [GE19]** (n = 2048, p = 10⁻³): the paper reports 20 million qubits and 8 hours
(5.1 h per run, 31 % retry risk). Our model gives 18 million and 6 hours. That is the same order
of magnitude, partly for different reasons: we use d = 31 instead of 27, 2d² · 1.5 instead of
2(d + 1)², 6,189 abstract qubits instead of a 14,238-patch layout, and factories at under 1 % of
the qubits instead of about 25 %. [G25] (preprint) lowers the count to under 1 million qubits in
under a week with newer arithmetic, yoked surface codes and magic-state cultivation.

**Newer estimates (2026 preprints, not peer reviewed), each with its own assumptions:**
- Pinnacle (Webster et al., arXiv:2602.11457): fewer than 100,000 qubits with qLDPC codes at
  p = 10⁻³, 1 µs cycle and 10 µs reaction time. qLDPC codes need longer-range connectivity than a
  2D nearest-neighbour grid.
- Mundada et al. (arXiv:2604.06319): 381,000 qubits and 9.2 days with a heterogeneous architecture
  (error rate and cycle time not stated in the abstract).
- Xue & Covey (arXiv:2605.03951): half a million qubits on a modular atomic processor with 1 ms
  measurement time.

How the two main papers reached their own numbers, step by step: [HOW_THE_PAPERS_COUNT.md](HOW_THE_PAPERS_COUNT.md).

## 5. Checking the rules with a base case

`verify_rules.py` builds or simulates small cases to test each rule the estimate uses. It runs in
about 35 s and uses fixed seeds, so a rerun gives the same numbers.

| Rule | How we test it | Result |
|---|---|---|
| R1. The cost is one block repeated, so it grows as a power of n | ccx count of real Qiskit circuits after two decomposition levels: Cuccaro adder (n = 4–32), cumulative (non-modular) multiplier (n = 4–16). Relative scaling only, not exact Toffoli counts | adder ∝ n^1.00, multiplier ∝ n^1.99, so exponentiation (2n multiplications) ∝ n^3 |
| R2. One surface-code patch uses 2d² − 1 physical qubits | Count the qubits in stim's rotated-memory circuits | 17, 49, 97, 161 for d = 3, 5, 7, 9: exact match |
| R3. Each step d → d + 2 divides the error by Λ ≈ 1 % / p | Surface-code memory, circuit-level depolarizing noise, pymatching, d = 3/5/7, up to 2 M shots | Λ = 9.94, 4.91, 2.46, 1.69, 1.13 at p = 0.1, 0.2, 0.4, 0.6, 1 %. The rule predicts 10, 5, 2.5, 1.7, 1.0 |
| R4. The risk of a run grows with its length | One d = 3 patch at p = 0.2 %, 3 to 81 rounds | P_fail = 0.30 %, 0.85 %, 2.5 %, 7.0 %: nearly a straight line |
| R5. Base case: RSA-2048 at p = 0.1 % | Fit the simulated d = 3, 5, 7 points and extend the line down to the 7.5·10⁻¹⁷ the run needs; Poisson bootstrap gives a band | d = 29 (68 % band 29–31); the [GE19] rule of thumb gives 31 |

What the check also shows:

- The rule of thumb is about **4× pessimistic** at d = 3, 5 and 7 and its slope matches the
  simulation, so it errs on the safe side by about one step of d.
- R5 extends three points over 13 orders of magnitude. At d = 7 and p = 0.1 % the simulation sees
  only 34 failures in 2 million shots, so Λ has a 68 % band of about 9.1–11.0.
- The simulation assumes ideal, independent noise. Real chips measure surface-code Λ = 1.4 to 2.14,
  which on our simulated Λ(p) curve corresponds to p ≈ 0.45 % (Λ = 2.14) to 0.8 % (Λ = 1.4); the formula Λ = 1 %/p gives 0.47 % to 0.71 %. Leakage,
  correlated errors and drift explain most of the gap.

## 6. Limitations

- Logical qubits, Toffolis and depth come from [GE19] formulas, not from a circuit we compiled.
- Runtime charges one 10 µs reaction time per layer of [GE19]'s measurement depth; it ignores
  factory throughput and the code-depth-limited lookups that [GE19] models.
- The error budget counts only the Q logical data qubits, not routing or factory patches.
- The factory term is crude: here it is under 1 % of the qubits; in [GE19] it is about a quarter.
- The rule of thumb is a fit; real devices have correlated errors, leakage and defects.

## References

- [B03] S. Beauregard, "Circuit for Shor's algorithm using 2n+3 qubits", Quantum Inf. Comput. 3, 175 (2003).
- [GE19] C. Gidney, M. Ekerå, "How to factor 2048 bit RSA integers in 8 hours using 20 million noisy qubits", Quantum 5, 433 (2021), arXiv:1905.09749.
- [G25] C. Gidney, "How to factor 2048 bit RSA integers with less than a million noisy qubits", arXiv:2505.15917 (2025, preprint).
- J. A. Smolin, G. Smith, A. Vargo, "Pretending to factor large numbers on a quantum computer", arXiv:1301.7007 (2013).
- Google Quantum AI, "Quantum error correction below the surface code threshold", Nature 638, 920 (2025), arXiv:2408.13687.
- T. He et al., "Experimental quantum error correction below the surface code threshold via all-microwave leakage suppression", Phys. Rev. Lett. 135, 260601 (2025).
- D. Bluvstein et al., "Architectural mechanisms of a universal fault-tolerant quantum computer", Nature (2025), doi:10.1038/s41586-025-09848-5, arXiv:2506.20661.
- Quantinuum, "Quantum error correction well below the color code threshold", arXiv:2610.24077 (2026, preprint).
- P. Webster et al., "The Pinnacle Architecture: Reducing the cost of breaking RSA-2048 to 100 000 physical qubits using quantum LDPC codes", arXiv:2602.11457 (2026, preprint).
- P. S. Mundada et al., "Heterogeneous architectures enable a 138x reduction in physical qubit requirements...", arXiv:2604.06319 (2026, preprint).
- T. Xue, J. P. Covey, "Factoring 2048 bit RSA integers with a half-million-qubit modular atomic processor", arXiv:2605.03951 (2026, preprint).

## License

MIT, see [LICENSE](LICENSE).
