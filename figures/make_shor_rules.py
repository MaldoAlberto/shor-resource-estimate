"""Figure "Do the rules hold? We checked": small simulations that confirm each rule behind the RSA-2048
estimates. Data: ../rules.json from verify_rules.py (Qiskit circuits for the
Toffoli counts; stim + pymatching surface-code memory with circuit-level noise, fixed seeds).
Hardware Lambda (surface code only): USTC PRL 2025, Google Nature 2025, Harvard Nature 2025.
Cream = intermediate category (post 2).
Usage (from the repo root): python3 figures/make_shor_rules.py <output folder>
"""
import json, pathlib, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

REPO = pathlib.Path(__file__).resolve().parents[1]          # repo root: data files live there
MASCOT = pathlib.Path.home() / "Downloads/QalyCat_character/qaly_happy.png"


def load_mascot():
    """The Qaly mascot is optional; the figure renders without it."""
    return np.asarray(Image.open(MASCOT)) if MASCOT.exists() else None


OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "out").resolve()
OUT.mkdir(parents=True, exist_ok=True)
R = json.loads((REPO / "rules.json").read_text())

INK, MUTED = "#13203A", "#5B6475"
TEAL, CORAL, GRAY, SAND = "#0A9396", "#C94A45", "#B7BDC8", "#E0A458"
BG = "#F6F4EF"  # cream = intermediate category
PCOL = {0.001: TEAL, 0.002: INK, 0.004: SAND, 0.01: CORAL}

plt.rcParams.update({"font.family": ["Helvetica Neue", "DejaVu Sans"], "font.size": 11, "text.color": INK,
                     "axes.edgecolor": MUTED, "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": BG,
                     "xtick.labelsize": 10, "ytick.labelsize": 10, "axes.labelsize": 10.5})


def render(footer, name):
    fig = plt.figure(figsize=(7.2, 9.0), dpi=300, facecolor=BG)
    fig.text(0.06, 0.955, "Do the rules hold?\nWe checked.", fontsize=22, fontweight="bold", va="top", linespacing=1.15)
    fig.text(0.06, 0.865, "Small simulations test each rule behind the RSA-2048 estimate", fontsize=12.5,
             color=MUTED, va="top")
    gs = fig.add_gridspec(3, 2, left=0.13, right=0.95, top=0.79, bottom=0.165, hspace=0.78, wspace=0.42)
    axes = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(6)]
    titles = ["1  Cost = one block, repeated", "2  Price of one square", "3  Bigger d, fewer errors",
              "4  Λ ≈ threshold / p", "5  Longer run, more risk", "6  Base case: RSA-2048"]
    for ax, t in zip(axes, titles):
        ax.set_title(t, loc="left", fontsize=12, fontweight="bold", color=INK, pad=8)

    # 1 Toffoli scaling of real Qiskit circuits
    ax = axes[0]
    rows = R["R1"]["rows"]
    n = np.array([r["n"] for r in rows])
    ax.loglog(n, [r["adder_toffolis"] for r in rows], "o-", color=TEAL, lw=1.8, ms=5, label="adder ~ n")
    m = [(r["n"], r["multiplier_toffolis"]) for r in rows if r["multiplier_toffolis"]]
    ax.loglog(*zip(*m), "s-", color=INK, lw=1.8, ms=5, label="multiplier ~ n²")
    ax.set_xticks([4, 8, 16, 32], ["4", "8", "16", "32"]); ax.minorticks_off()
    ax.set_xlabel("bits n"); ax.set_ylabel("Toffoli gates")
    ax.text(4, 1.6e3, "multiplier ~ n²", fontsize=9, color=INK)
    ax.text(32, 120, "adder ~ n", fontsize=9, color=TEAL, ha="right")
    ax.set_ylim(4, 2e4)

    # 2 qubits per patch
    ax = axes[1]
    dd = np.linspace(3, 9, 50)
    ax.plot(dd, 2 * dd ** 2 - 1, color=GRAY, lw=1.6, label="rule 2d² − 1")
    ax.plot([r["d"] for r in R["R2"]], [r["qubits"] for r in R["R2"]], "o", color=TEAL, ms=7, label="counted in stim")
    for r in R["R2"]:
        ax.text(r["d"] + (0.3 if r["d"] == 3 else -0.25), r["qubits"] + 9, str(r["qubits"]), fontsize=9.5,
                ha="left" if r["d"] == 3 else "right")
    ax.set_xticks([3, 5, 7, 9]); ax.set_ylim(0, 200)
    ax.set_xlabel("code distance d"); ax.set_ylabel("physical qubits")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    # 3 logical error vs d for several p
    ax = axes[2]
    for row in R["R3"]:
        p = row["p"]
        if p not in PCOL:
            continue
        ds = [int(d) for d in row["pL"]]
        v = [row["pL"][str(d)]["per_round"] for d in ds]
        ax.semilogy(ds, v, "o-", color=PCOL[p], lw=1.8, ms=5)
        ax.text(7.3, v[-1], f"p = {p * 100:g} %", fontsize=8.5, color=PCOL[p], va="center")
    ax.set_xlim(2.5, 9.6); ax.set_xticks([3, 5, 7])
    ax.set_xlabel("code distance d"); ax.set_ylabel("logical error / round")

    # 4 Lambda vs p, rule and hardware
    ax = axes[3]
    ps = np.array([r["p"] for r in R["R3"]]); lam = np.array([r["Lambda"] for r in R["R3"]])
    pp = np.logspace(np.log10(8e-4), np.log10(1.2e-2), 50)
    ax.loglog(pp, 0.01 / pp, color=GRAY, lw=1.6, label="rule: 1 % / p")
    ax.loglog(ps, lam, "o", color=TEAL, ms=7, label="simulated")
    for lab, v, c in [("Google, Harvard 2.14", 2.14, INK), ("USTC 1.4", 1.4, INK)]:
        ax.axhline(v, color=c, lw=0.8, ls=":")
        ax.text(8.6e-4, v * 1.08, lab, fontsize=8, ha="left", va="bottom", color=c)
    ax.set_ylim(0.8, 20); ax.set_yticks([1, 2, 5, 10], ["1", "2", "5", "10"]); ax.minorticks_off()
    ax.set_xticks([1e-3, 1e-2], ["0.1 %", "1 %"])
    ax.set_xlabel("gate error p"); ax.set_ylabel("Λ (÷ per d + 2)")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right", bbox_to_anchor=(1.0, 0.88))

    # 5 failure vs rounds
    ax = axes[4]
    r4 = R["R4"]
    x = np.array([r["rounds"] for r in r4]); y = np.array([r["p_fail"] for r in r4])
    k = np.sum(x * y) / np.sum(x * x)
    ax.plot([0, 90], [0, 90 * k], color=GRAY, lw=1.6, label="straight line")
    ax.plot(x, y, "o", color=TEAL, ms=7, label="simulated, d = 3")
    ax.set_xlabel("rounds"); ax.set_ylabel("chance the run fails")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    # 6 base case: extrapolate measured curve to RSA-2048
    ax = axes[5]
    r5 = R["R5"]
    dd = np.arange(3, 36, 2)
    meas = next(r for r in R["R3"] if r["p"] == 0.001)
    ax.semilogy(dd, r5["pL_d3"] * r5["Lambda"] ** (-(dd - 3) / 2), "--", color=TEAL, lw=1.5, label="simulated, extended")
    ax.semilogy(dd, 0.1 * 0.1 ** ((dd + 1) / 2), ":", color=CORAL, lw=1.6, label="formula")
    ax.semilogy([3, 5, 7], [meas["pL"][k]["per_round"] for k in ("3", "5", "7")], "o", color=TEAL, ms=6)
    ax.axhline(r5["need"], color=MUTED, lw=1, ls="-")
    ax.text(3, r5["need"] * 4, "RSA-2048 needs", fontsize=8.5, color=MUTED, va="bottom")
    ax.plot([r5["d_measured_curve"]], [r5["need"]], "o", color=INK, ms=7)
    ax.text(r5["d_measured_curve"] - 1.5, r5["need"] * 0.1,
            f"d = {r5['d_measured_curve']}\n(band {r5['d_band'][0]}–{r5['d_band'][1]};\nformula {r5['d_formula']})",
            fontsize=9, ha="right", va="top", linespacing=1.2)
    ax.set_ylim(1e-25, 1e-2); ax.set_yticks([1e-3, 1e-10, 1e-17])
    ax.set_xlabel("code distance d  (p = 0.1 %)"); ax.set_ylabel("logical error / round")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    fig.text(0.06, 0.048,
             "Panel 1: real Qiskit adders and multipliers. Panels 2–6: surface-code memory, circuit-level noise\n"
             "(stim + pymatching), up to 2 million shots, fixed seeds. Panel 1 shows relative scaling only. The formula\n"
             "is ~4× pessimistic with the same slope: simulation d = 29 (68 % band 29–31), formula d = 31.\n"
             "Hardware lines: surface-code Λ, USTC PRL 2025, Google Nature 2025, Harvard Nature 2025.",
             fontsize=7.8, color=MUTED, linespacing=1.4)
    cat = load_mascot()
    w = 0.13
    cax = fig.add_axes([0.82, 0.855, w, w * (396 / 439) * (7.2 / 9.0)])
    cax.axis("off")
    if cat is not None:
        cax.imshow(cat)
    fig.text(0.80, 0.915, "Don't trust a rule\nyou haven't run!", fontsize=9.5, ha="right", va="center",
             color=INK, linespacing=1.3,
             bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.9", fc="white", ec=TEAL, lw=1.5))
    fig.text(0.06, 0.015, footer, fontsize=11, color=MUTED)
    fig.text(0.95, 0.015, "@qalycat", fontsize=11, color=MUTED, ha="right")
    png = OUT / f"{name}.png"
    fig.savefig(png, facecolor=BG); plt.close(fig)
    Image.open(png).convert("RGB").resize((1080, 1350), Image.LANCZOS).save(OUT / f"{name}.jpg", quality=92)
    return png


print(render("Code: link in bio", "shor-rules-en-instagram"))
print(render("Code: github.com/MaldoAlberto/shor-resource-estimate", "shor-rules-en-linkedin"))
