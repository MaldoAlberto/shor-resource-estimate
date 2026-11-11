"""Shor's algorithm for N = 15 built from real modular arithmetic (Beauregard 2003), not a compiled oracle.

Why: small Shor demos often replace the modular multiplication by a precomputed permutation (which
needs the full classical multiplication table) or simplify it with knowledge of the period (the
problem described by Smolin, Smith, Vargo, arXiv:1301.7007). Here every controlled multiplication is
built from Draper QFT adders -> modular adder -> modular multiplier -> controlled U_a, for ANY a and N,
and we keep all of them, even when a^(2^j) mod N = 1.

Checks:
  1. U_a is verified on every input x < N, with control = 1 (x -> a x mod N, work qubits back to 0)
     and control = 0 (nothing changes).
  2. The full circuit is simulated; the counting register shows the period; 15 = 3 x 5.
  3. Gate counts after compiling to {cx, rz, sx, x} (optimization level 1, same as shor_estimate.py).

Usage: python3 shor_arithmetic.py [--count 8] [--json arith.json]
"""
import argparse, json, math, time

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit.circuit.library import QFTGate

from shor_estimate import check_inputs, factors_from_peaks
from qiskit.quantum_info import Statevector


def phi_add(nb, a):
    """Add the constant a to a register already in the Fourier basis (Draper): one phase per qubit."""
    qc = QuantumCircuit(nb, name=f"φADD({a})")
    for j in range(nb):
        qc.p(2 * math.pi * a * 2 ** j / 2 ** nb, j)
    return qc


def c_gate(circ, k):
    """Circuit -> gate with k extra control qubits (controls come first)."""
    return circ.to_gate().control(k)


def phi_add_mod(n, a, N):
    """Beauregard Fig. 5: doubly-controlled  |b> -> |(a + b) mod N>  in Fourier space.
    Qubits: [c1, c2, b_0..b_n, anc]  (b has n + 1 qubits, b < N)."""
    nb = n + 1
    c = QuantumRegister(2, "c"); b = QuantumRegister(nb, "b"); anc = QuantumRegister(1, "anc")
    qc = QuantumCircuit(c, b, anc, name=f"φADDMOD({a})")
    qft, iqft = QFTGate(nb), QFTGate(nb).inverse()
    msb = b[nb - 1]
    qc.append(c_gate(phi_add(nb, a), 2), [c[0], c[1], *b])
    qc.append(phi_add(nb, -N).to_gate(), b)
    qc.append(iqft, b); qc.cx(msb, anc[0]); qc.append(qft, b)
    qc.append(c_gate(phi_add(nb, N), 1), [anc[0], *b])
    qc.append(c_gate(phi_add(nb, -a), 2), [c[0], c[1], *b])
    qc.append(iqft, b); qc.x(msb); qc.cx(msb, anc[0]); qc.x(msb); qc.append(qft, b)
    qc.append(c_gate(phi_add(nb, a), 2), [c[0], c[1], *b])
    return qc


def c_mult_mod(n, a, N):
    """Controlled  |x>|b> -> |x>|(b + a x) mod N>.  Qubits: [ctrl, x_0..x_{n-1}, b_0..b_n, anc]."""
    nb = n + 1
    ctrl = QuantumRegister(1, "ctrl"); x = QuantumRegister(n, "x"); b = QuantumRegister(nb, "b")
    anc = QuantumRegister(1, "anc")
    qc = QuantumCircuit(ctrl, x, b, anc, name=f"CMULT({a})")
    qc.append(QFTGate(nb), b)
    for i in range(n):
        qc.append(phi_add_mod(n, (a * 2 ** i) % N, N).to_gate(), [ctrl[0], x[i], *b, anc[0]])
    qc.append(QFTGate(nb).inverse(), b)
    return qc


def c_u(n, a, N):
    """Controlled U_a: |x> -> |a x mod N>, work register b and anc restored to |0>.
    Qubits: [ctrl, x_0..x_{n-1}, b_0..b_n, anc]."""
    nb = n + 1
    ctrl = QuantumRegister(1, "ctrl"); x = QuantumRegister(n, "x"); b = QuantumRegister(nb, "b")
    anc = QuantumRegister(1, "anc")
    qc = QuantumCircuit(ctrl, x, b, anc, name=f"C-U({a})")
    allq = [ctrl[0], *x, *b, anc[0]]
    qc.append(c_mult_mod(n, a, N).to_gate(), allq)
    for i in range(n):
        qc.cswap(ctrl[0], x[i], b[i])
    qc.append(c_mult_mod(n, pow(a, -1, N), N).inverse().to_gate(), allq)
    return qc


def check_cu(n, a, N):
    """Run C-U_a on every basis input x < N with control = 1 and control = 0; return the failures."""
    circ = c_u(n, a, N)
    bad = []
    for ctrl in (1, 0):
        for xval in range(N):
            idx = ctrl | (xval << 1)                          # ctrl is qubit 0, x in qubits 1..n
            sv = Statevector.from_int(idx, 2 ** circ.num_qubits).evolve(circ)
            out = int(np.argmax(np.abs(sv.data)))
            want = ctrl | ((((a * xval) % N) if ctrl else xval) << 1)
            if out != want or abs(sv.data[out] - 1) > 1e-6:  # also catches a wrong phase
                bad.append((ctrl, xval, out, want))
    return bad


def shor_circuit(N, a, n_count):
    n = N.bit_length()
    cnt = QuantumRegister(n_count, "count"); x = QuantumRegister(n, "x")
    b = QuantumRegister(n + 1, "b"); anc = QuantumRegister(1, "anc")
    qc = QuantumCircuit(cnt, x, b, anc)
    qc.h(cnt)
    qc.x(x[0])                                                # x = 1
    for j in range(n_count):                                  # no shortcut, even when a^(2^j) mod N = 1
        qc.append(c_u(n, pow(a, 2 ** j, N), N).to_gate(), [cnt[j], *x, *b, anc[0]])
    qc.append(QFTGate(n_count).inverse(), cnt)
    return qc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=15)
    ap.add_argument("--a", type=int, default=7)
    ap.add_argument("--count", type=int, default=8)
    ap.add_argument("--json")
    args = ap.parse_args()
    N, a, n = args.N, args.a, args.N.bit_length()
    check_inputs(N, a)
    t0 = time.time()

    bad = check_cu(n, a, N)
    print(f"1. C-U_{a} checked on all {N} inputs x < {N}, control = 1 and 0: {'OK' if not bad else f'FAILED {bad}'}")

    qc = shor_circuit(N, a, args.count)
    probs = Statevector(qc).probabilities(list(range(args.count)))
    peaks = {int(i): float(p) for i, p in enumerate(probs) if p > 1e-3}
    res = factors_from_peaks(peaks, args.count, a, N)
    if res is None:
        raise SystemExit(f"no usable period for a = {a} (odd r or a^(r/2) = -1 mod N); try another a")
    r, f1, f2 = res
    print(f"2. {qc.num_qubits} qubits simulated. peaks: "
          + ", ".join(f"{y}: {p:.3f}" for y, p in sorted(peaks.items())) + f"   r = {r}  ->  {N} = {f1} x {f2}")

    tqc = transpile(qc, basis_gates=["cx", "rz", "sx", "x"], optimization_level=1, seed_transpiler=7)
    ops = dict(tqc.count_ops())
    counts = {"qubits": qc.num_qubits, "cx": ops.get("cx", 0),
              "one_qubit": sum(v for k, v in ops.items() if k in ("rz", "sx", "x")), "depth": tqc.depth()}
    print(f"3. compiled to cx/rz/sx/x: {counts['qubits']} qubits, {counts['cx']:,} CX, "
          f"{counts['one_qubit']:,} one-qubit gates, depth {counts['depth']:,}")
    one = transpile(c_u(n, a, N), basis_gates=["cx", "rz", "sx", "x"], optimization_level=1, seed_transpiler=7)
    counts["cx_one_multiplication"] = dict(one.count_ops()).get("cx", 0)
    print(f"   one controlled x{a} mod {N}: {counts['cx_one_multiplication']:,} CX "
          f"(all {args.count} kept, including the ones equal to x1)")
    print(f"   ({time.time() - t0:.0f} s)")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"N": N, "a": a, "n_count": args.count, "cu_check_failures": bad, "peaks": peaks,
                       "period": r, "factors": [f1, f2], "counts": counts}, fh, indent=1)


if __name__ == "__main__":
    main()
