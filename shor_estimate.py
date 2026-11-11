"""How many qubits and gates does Shor's algorithm need?

Two parts, same method as the DQC tutorial (known answer first, then scale up):

  Part A  Reference circuit for N = 15, a = 7 with a PRECOMPUTED-PERMUTATION ORACLE: each controlled
          multiplication is a permutation matrix built from the full classical multiplication table
          (exponential classical cost), and the compiler drops the multiplications that equal 1.
          It is kept only as a comparison point for shor_arithmetic.py, which builds the same circuit
          from real modular arithmetic.

  Part B  Scale up: estimate logical qubits, Toffoli gates, surface-code distance, physical qubits
          and runtime for an n-bit modulus (up to RSA-2048), from published formulas plus a simple,
          explicit surface-code model. These are ESTIMATES; the assumptions are printed.

Usage:
  python3 shor_estimate.py                     # Part A + table for several sizes
  python3 shor_estimate.py --skip-a --bits 2048 --p 1e-3
  python3 shor_estimate.py --json results.json # save everything (used by the figure scripts)

Sources for Part B:
  [B03] S. Beauregard, "Circuit for Shor's algorithm using 2n+3 qubits", QIC 3, 175 (2003).
  [GE19] C. Gidney, M. Ekera, "How to factor 2048 bit RSA integers in 8 hours using 20 million
         noisy qubits", Quantum 5, 433 (2021), arXiv:1905.09749. Logical-error rule of thumb
         0.1 (p / 1 %)^((d+1)/2) and 10 us reaction-limited layers are taken from this paper.
  [G25] C. Gidney, "How to factor 2048 bit RSA integers with less than a million noisy qubits",
        arXiv:2505.15917 (2025, preprint).
"""
import argparse, json, math
from fractions import Fraction

import numpy as np

# ===================================================================== shared classical post-processing


def factors_from_peaks(peaks, n_count, a, N):
    """Continued fractions on the measured peaks -> period r -> factors.
    Returns (r, f1, f2), or None if no peak gives a usable even period (then rerun with another a)."""
    for y in sorted(peaks, key=peaks.get, reverse=True):
        den = Fraction(y, 2 ** n_count).limit_denominator(N).denominator
        for r in (den * k for k in range(1, 5)):              # the fraction may give a divisor of r
            if pow(a, r, N) != 1:
                continue
            if r % 2 or pow(a, r // 2, N) == N - 1:            # unusable period: try another a
                break
            f1, f2 = math.gcd(pow(a, r // 2) - 1, N), math.gcd(pow(a, r // 2) + 1, N)
            if 1 < f1 < N:
                return r, f1, f2
            break
    return None


def check_inputs(N, a):
    if N % 2 == 0 or N < 15:
        raise ValueError("N must be odd and at least 15")
    g = math.gcd(a, N)
    if g != 1:
        raise ValueError(f"gcd(a, N) = {g}: a already shares a factor with N, no quantum step needed")


# ===================================================================== Part A: precomputed oracle


def part_a(N=15, a=7, n_count=8):
    from qiskit import QuantumCircuit, transpile
    from qiskit.circuit.library import QFTGate, UnitaryGate
    from qiskit.quantum_info import Statevector

    check_inputs(N, a)
    n_work = N.bit_length()

    def mult_gate(k):
        """|y> -> |a^k * y mod N> for y < N (identity for y >= N): a permutation built from the full table."""
        m = pow(a, k, N)
        dim = 2 ** n_work
        U = np.zeros((dim, dim))
        for y in range(dim):
            U[(m * y) % N if y < N else y, y] = 1
        return UnitaryGate(U, label=f"x{m} mod {N}").control(1)

    qc = QuantumCircuit(n_count + n_work)
    qc.h(range(n_count))
    qc.x(n_count)                                             # work register starts in |1>
    for j in range(n_count):                                  # controlled a^(2^j) mod N
        qc.append(mult_gate(2 ** j), [j] + list(range(n_count, n_count + n_work)))
    qc.append(QFTGate(n_count).inverse(), range(n_count))

    # same compiler settings as shor_arithmetic.py
    tqc = transpile(qc, basis_gates=["cx", "rz", "sx", "x"], optimization_level=1, seed_transpiler=7)
    ops = dict(tqc.count_ops())
    counts = {"qubits": qc.num_qubits, "cx": ops.get("cx", 0),
              "one_qubit": sum(v for k, v in ops.items() if k in ("rz", "sx", "x")), "depth": tqc.depth()}
    # one controlled multiplication by 7 alone, for a like-for-like comparison with shor_arithmetic.py
    one = QuantumCircuit(1 + n_work); one.append(mult_gate(1), range(1 + n_work))
    cx_one = dict(transpile(one, basis_gates=["cx", "rz", "sx", "x"], optimization_level=1,
                            seed_transpiler=7).count_ops()).get("cx", 0)

    probs = Statevector(qc).probabilities(list(range(n_count)))
    peaks = {int(i): float(p) for i, p in enumerate(probs) if p > 1e-6}
    res = factors_from_peaks(peaks, n_count, a, N)
    if res is None:
        raise RuntimeError(f"no usable period for a = {a}; try another a")
    r, f1, f2 = res
    return {"N": N, "a": a, "n_count": n_count, "counts": counts, "cx_one_multiplication": cx_one,
            "peaks": peaks, "period": r, "factors": [f1, f2]}


# ===================================================================== Part B: scale up, estimated

CYCLE_US = 1.0          # one surface-code round = 1 microsecond (superconducting-like, as in GE19)
REACTION_ROUNDS = 10    # one logical layer waits one reaction time = 10 rounds (GE19: 10 us at 1 us)
P_TH = 1e-2             # threshold in the GE19 rule of thumb (a property of code + noise + decoder)
FAIL_BUDGET = 0.01      # OUR choice: a 1 % chance that the whole run fails (GE19 accepts 31 % and retries)
ROUTING = 1.5           # OUR assumption: extra patches for routing / lattice surgery


def logical_qubits_beauregard(n):            # [B03], minimal-qubit textbook circuit
    return 2 * n + 3


def logical_qubits_ge(n):                    # [GE19] abstract formula
    return 3 * n + 0.002 * n * math.log2(n)


def toffolis_ge(n):                          # [GE19] Toffoli (+T/2) count, windowed arithmetic
    return 0.3 * n ** 3 + 0.0005 * n ** 3 * math.log2(n)


def meas_depth_ge(n):                        # [GE19] measurement depth (reaction-limited layers)
    return 500 * n ** 2 + n ** 2 * math.log2(n)


def p_logical(d, p):
    """Logical error per logical qubit per round, GE19 rule of thumb: 0.1 (p / p_th)^((d+1)/2)."""
    return 0.1 * (p / P_TH) ** ((d + 1) / 2)


def estimate(n, p=1e-3, cycle_us=CYCLE_US):
    q_log = logical_qubits_ge(n)
    tof = toffolis_ge(n)
    depth = meas_depth_ge(n)
    rounds = depth * REACTION_ROUNDS          # every layer waits one reaction time
    for d in range(3, 202, 2):                # smallest odd distance that meets the budget
        if q_log * rounds * p_logical(d, p) < FAIL_BUDGET:
            break
    else:                                     # p too close to threshold: no distance <= 201 works
        return None
    phys_data = q_log * 2 * d ** 2 * ROUTING  # 2d^2 - 1 ~ 2d^2 physical qubits per patch, x routing
    phys_fact = 4 * 12 * 2 * d ** 2           # OUR crude guess: 4 factories of 12 patches each
    return {"n": n, "p": p, "logical_qubits": q_log, "logical_qubits_min": logical_qubits_beauregard(n),
            "toffolis": tof, "meas_depth": depth, "distance": d,
            "physical_qubits": phys_data + phys_fact, "hours": rounds * cycle_us / 3.6e9}


PUBLISHED = [  # n = 2048, p = 1e-3
    {"label": "Gidney-Ekera 2019", "n": 2048, "physical_qubits": 20e6, "hours": 8},
    {"label": "Gidney 2025 (preprint)", "n": 2048, "physical_qubits": 1e6, "hours": 7 * 24, "upper_bound": True},
]


def fmt(x):
    return f"{x:,.0f}" if x < 1e6 else f"{x:.2e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="*", default=[8, 64, 256, 1024, 2048, 4096])
    ap.add_argument("--p", type=float, default=1e-3, help="physical error rate per gate")
    ap.add_argument("--json", help="save results to this file")
    ap.add_argument("--skip-a", action="store_true")
    args = ap.parse_args()

    out = {}
    if not args.skip_a:
        A = part_a()
        out["part_a"] = A
        c = A["counts"]
        print("Part A  N = 15, a = 7 with a precomputed-permutation oracle (comparison only)")
        print(f"  qubits: {c['qubits']}  CX: {c['cx']}  1-qubit gates: {c['one_qubit']}  depth: {c['depth']}"
              f"  (one controlled x7 mod 15: {A['cx_one_multiplication']} CX)")
        top = sorted(A["peaks"].items(), key=lambda kv: -kv[1])
        print("  peaks (y: prob):", ", ".join(f"{y}: {p:.3f}" for y, p in top))
        print(f"  period r = {A['period']}  ->  15 = {A['factors'][0]} x {A['factors'][1]}\n")

    print(f"Part B  estimates, p = {args.p:g}, round = {CYCLE_US} us, layer = {REACTION_ROUNDS} rounds "
          f"(reaction-limited), failure budget = {FAIL_BUDGET:.0%}")
    print(f"  {'bits':>5} {'log.q(min)':>10} {'log.q':>8} {'Toffoli':>10} {'d':>3} {'phys.q':>10} {'hours':>8}")
    rows = []
    for n in args.bits:
        e = estimate(n, args.p)
        if e is None:
            print(f"  {n:>5}  no code distance <= 201 meets the budget at p = {args.p:g} (too close to threshold)")
            continue
        rows.append(e)
        print(f"  {n:>5} {e['logical_qubits_min']:>10,} {e['logical_qubits']:>8,.0f} {e['toffolis']:>10.2e} "
              f"{e['distance']:>3} {fmt(e['physical_qubits']):>10} {e['hours']:>8.2f}")
    print("  published (n = 2048, p = 1e-3):",
          "; ".join(f"{r['label']}: {'<' if r.get('upper_bound') else ''}{r['physical_qubits']:.0e} qubits, "
                    f"{'<' if r.get('upper_bound') else ''}{r['hours']:.0f} h" for r in PUBLISHED))
    out["part_b"] = {"assumptions": {"cycle_us": CYCLE_US, "reaction_rounds": REACTION_ROUNDS, "p_th": P_TH,
                                     "fail_budget": FAIL_BUDGET, "routing": ROUTING,
                                     "factories": "4 x 12 patches (crude)"},
                     "rows": rows, "published": PUBLISHED}
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
        print(f"\nsaved {args.json}")


if __name__ == "__main__":
    main()
