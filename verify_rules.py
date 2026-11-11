"""Check, by building and simulating small cases, the five rules every RSA-2048 estimate relies on.

  R1  The algorithm is one block repeated many times, so its cost scales with a power of n.
      We build real adders and multipliers in Qiskit and count their ccx gates after two
      decomposition levels. These are RELATIVE SCALING counts, not exact Toffoli counts (the HRS
      multiplier is not modular): adder ~ n, multiplier ~ n^2, so exponentiation (~2n multiplications) ~ n^3.
  R2  One rotated surface-code logical qubit uses 2d^2 - 1 physical qubits.
      We count the qubits in stim's generated memory circuits.
  R3  Below threshold, every step d -> d+2 divides the logical error by Lambda ~ p_th / p.
      We simulate surface-code memory with circuit-level noise (stim) and decode it (pymatching).
  R4  The risk of a long run grows with rounds: P_fail(R rounds) ~ R * p_L (while small).
      We simulate one patch for more and more rounds.
  R5  Base case: extend the simulated curve (not the formula) to the error RSA-2048 needs, with a
      Poisson-bootstrap band on the fitted slope, and read off the code distance.

Usage (run from the repo folder, it imports shor_estimate.py):
  python3 verify_rules.py [--quick] [--json rules.json]
Needs: see requirements.txt (stim==1.14.0 and pymatching==2.2.1 are the last wheels for Python 3.9).
"""
import argparse, json, time

import numpy as np

# ===================================================================== R1: cost of repeated blocks


def r1_toffoli_scaling(ns=(4, 8, 16, 32)):
    from qiskit.circuit.library import CDKMRippleCarryAdder, HRSCumulativeMultiplier

    def tof(circ):
        # two decomposition levels expose the Toffolis (ccx) inside the MAJ/UMA blocks of every adder
        return dict(circ.decompose(reps=2).count_ops()).get("ccx", 0)

    rows = []
    for n in ns:
        add = tof(CDKMRippleCarryAdder(n, kind="half"))
        mul = tof(HRSCumulativeMultiplier(n, adder=CDKMRippleCarryAdder(n, kind="half"))) if n <= 16 else None
        rows.append({"n": n, "adder_toffolis": add, "multiplier_toffolis": mul})
    nn = np.array([r["n"] for r in rows], float)
    a = np.array([r["adder_toffolis"] for r in rows], float)
    k_add = np.polyfit(np.log(nn), np.log(a), 1)[0]
    mrows = [r for r in rows if r["multiplier_toffolis"]]
    k_mul = np.polyfit(np.log([r["n"] for r in mrows]), np.log([r["multiplier_toffolis"] for r in mrows]), 1)[0]
    return {"rows": rows, "exponent_adder": k_add, "exponent_multiplier": k_mul}


# ===================================================================== R2-R4: surface code with stim


def memory_circuit(d, p, rounds):
    import stim
    return stim.Circuit.generated("surface_code:rotated_memory_z", distance=d, rounds=rounds,
                                  after_clifford_depolarization=p, before_round_data_depolarization=p,
                                  before_measure_flip_probability=p, after_reset_flip_probability=p)


def logical_failure(d, p, rounds, max_shots, max_errors=200, batch=100_000):
    """Monte Carlo: fraction of shots where the decoded logical observable is wrong."""
    import pymatching
    c = memory_circuit(d, p, rounds)
    sampler = c.compile_detector_sampler(seed=7 + 1000 * d + rounds)   # fixed seed: rerun gives same numbers
    m = pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
    shots = errs = 0
    while shots < max_shots and errs < max_errors:
        b = min(batch, max_shots - shots)
        det, obs = sampler.sample(b, separate_observables=True)
        errs += int(np.sum(m.decode_batch(det)[:, 0] != obs[:, 0]))
        shots += b
    return errs, shots


def per_round(p_shot, rounds):
    return 1 - (1 - p_shot) ** (1 / rounds)


def r2_qubits(ds=(3, 5, 7, 9)):
    # num_qubits counts unused grid indices too, so count the qubits that have coordinates (= are used)
    return [{"d": d, "qubits": len(memory_circuit(d, 1e-3, d).get_final_qubit_coordinates()), "rule": 2 * d * d - 1}
            for d in ds]


def r3_lambda(ps, ds, max_shots):
    out = []
    for p in ps:
        row = {"p": p, "pL": {}}
        for d in ds:
            e, s = logical_failure(d, p, rounds=d, max_shots=max_shots)
            row["pL"][d] = {"errors": e, "shots": s, "per_round": per_round(e / s, d) if e else None}
        pl = [(d, row["pL"][d]["per_round"]) for d in ds if row["pL"][d]["per_round"]]
        if len(pl) >= 2:
            slope = np.polyfit([d for d, _ in pl], np.log([v for _, v in pl]), 1)[0]
            row["Lambda"] = float(np.exp(-2 * slope))            # factor per d -> d+2
        out.append(row)
    return out


def r4_rounds(d, p, rounds_list, max_shots):
    out = []
    for r in rounds_list:
        e, s = logical_failure(d, p, rounds=r, max_shots=max_shots, max_errors=2000)
        out.append({"rounds": r, "errors": e, "shots": s, "p_fail": e / s})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="fewer shots, for a fast check")
    ap.add_argument("--json")
    a = ap.parse_args()
    shots = 200_000 if a.quick else 2_000_000
    out = {}

    t = time.time()
    R1 = r1_toffoli_scaling()
    out["R1"] = R1
    print("R1  ccx count of real Qiskit circuits (Cuccaro adder, cumulative multiplier); relative scaling only")
    for r in R1["rows"]:
        print(f"    n = {r['n']:>3}: adder {r['adder_toffolis']:>5}   multiplier {r['multiplier_toffolis'] or '-':>6}")
    print(f"    fitted exponent: adder n^{R1['exponent_adder']:.2f}, multiplier n^{R1['exponent_multiplier']:.2f}"
          f"  -> exponentiation (2n multiplications) ~ n^{R1['exponent_multiplier'] + 1:.2f}")

    R2 = r2_qubits()
    out["R2"] = R2
    print("\nR2  physical qubits per rotated surface-code patch (stim)")
    for r in R2:
        print(f"    d = {r['d']}: {r['qubits']} qubits, rule 2d^2-1 = {r['rule']}  {'OK' if r['qubits'] == r['rule'] else 'MISMATCH'}")

    R3 = r3_lambda(ps=[1e-3, 2e-3, 4e-3, 6e-3, 1e-2], ds=[3, 5, 7], max_shots=shots)
    out["R3"] = R3
    print("\nR3  logical error per round (circuit-level noise, pymatching), and Lambda per d -> d+2")
    for row in R3:
        cells = "  ".join(f"d={d}: {v['per_round']:.1e} ({v['errors']}/{v['shots']})" if v["per_round"]
                          else f"d={d}: <{1 / v['shots']:.0e} (0/{v['shots']})" for d, v in row["pL"].items())
        lam = row.get("Lambda")
        print(f"    p = {row['p']:.0e}: {cells}   Lambda = {lam:.2f}  (rule p_th/p with p_th = 1 %: {0.01 / row['p']:.1f})"
              if lam else f"    p = {row['p']:.0e}: {cells}")
        for d, v in row["pL"].items():
            if v["per_round"]:
                v["formula"] = 0.1 * (row["p"] / 0.01) ** ((d + 1) / 2)

    R4 = r4_rounds(d=3, p=2e-3, rounds_list=[3, 9, 27, 81], max_shots=shots)
    out["R4"] = R4
    print("\nR4  one d = 3 patch at p = 2e-3: failure vs number of rounds")
    for r in R4:
        print(f"    {r['rounds']:>3} rounds: P_fail = {r['p_fail']:.2e}  ({r['errors']}/{r['shots']})")
    # R5: worked base case. Pick d for RSA-2048 with the MEASURED curve at p = 1e-3, not the formula.
    import shor_estimate as se
    row = next(r for r in R3 if r["p"] == 1e-3)
    pl3, lam = row["pL"][3]["per_round"], row["Lambda"]
    E = se.estimate(2048, 1e-3)
    # budget per logical qubit per round; does not depend on d (rounds = layers x reaction time)
    need = se.FAIL_BUDGET / (E["logical_qubits"] * E["meas_depth"] * se.REACTION_ROUNDS)
    ds = np.array([3, 5, 7])

    def d_needed(errs):
        pl = [per_round(e / row["pL"][d]["shots"], d) for d, e in zip(ds, errs)]
        slope, icpt = np.polyfit(ds, np.log(pl), 1)
        d = 3
        while np.exp(icpt + slope * d) > need and d < 99:
            d += 2
        return d, float(np.exp(-2 * slope)), float(np.exp(icpt + slope * 3))

    errs = np.array([row["pL"][d]["errors"] for d in ds])
    d_c, lam_c, pl3_c = d_needed(errs)
    rng = np.random.default_rng(7)
    boot = [d_needed(rng.poisson(errs))[0] for _ in range(2000)]
    d_lo, d_hi = (int(v) for v in np.percentile(boot, [16, 84]))
    lam_boot = [d_needed(rng.poisson(errs))[1] for _ in range(500)]
    lam_lo, lam_hi = np.percentile(lam_boot, [16, 84])
    out["R5"] = {"pL_d3": pl3_c, "Lambda": lam_c, "Lambda_band": [float(lam_lo), float(lam_hi)], "need": need,
                 "d_measured_curve": d_c, "d_band": [d_lo, d_hi], "d_formula": E["distance"]}
    print(f"\nR5  base case: RSA-2048 at p = 1e-3 needs p_L < {need:.1e} per qubit per round")
    print(f"    fitted simulated curve: Lambda = {lam_c:.2f} (68 % band {lam_lo:.1f}-{lam_hi:.1f})  ->  d = {d_c}"
          f" (68 % band {d_lo}-{d_hi});  GE19 rule of thumb gives d = {E['distance']}")
    for dd in (3, 5, 7):
        v = row["pL"][dd]
        print(f"    d = {dd}: simulated {v['per_round']:.1e}, formula {v['formula']:.1e} "
              f"(formula / simulation = {v['formula'] / v['per_round']:.1f})")

    print(f"\n    done in {time.time() - t:.0f} s")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1, default=float)
        print(f"    saved {a.json}")


if __name__ == "__main__":
    main()
