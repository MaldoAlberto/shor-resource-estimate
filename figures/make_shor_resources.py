"""Figure "How big is Shor's algorithm?": factor 15 for real (measured gate counts + simulated
peaks), then estimate logical qubits, Toffolis, code distance and physical qubits up to RSA-2048.
Data: ../shor_estimate.py (Part A measured, Part B estimated + published).
Cream background = intermediate category. English, LinkedIn + Instagram (Qaly, "link in bio").
Usage (from the repo root): python3 figures/make_shor_resources.py <output folder>     (needs matplotlib, numpy, Pillow, qiskit)
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


sys.path.insert(0, str(REPO))
import shor_estimate as se  # noqa: E402

OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "out").resolve()
OUT.mkdir(parents=True, exist_ok=True)

INK, MUTED = "#13203A", "#5B6475"
TEAL, CORAL, GRAY = "#0A9396", "#C94A45", "#B7BDC8"
BG = "#F6F4EF"  # cream = intermediate category (post 2 of the Shor series)

ARITH = REPO / "arith.json"       # from shor_arithmetic.py (real arithmetic)
A = json.loads(ARITH.read_text()); A["peaks"] = {int(k): v for k, v in A["peaks"].items()}
O = se.part_a()                                                          # precomputed-permutation oracle, for comparison
O["counts"]["cx_one_multiplication"] = O["cx_one_multiplication"]
bits = np.unique(np.round(np.logspace(np.log10(8), np.log10(4096), 40)).astype(int))
E3 = [se.estimate(int(n), 1e-3) for n in bits]
E4 = [se.estimate(int(n), 1e-4) for n in bits]
col = lambda E, k: np.array([e[k] for e in E])

plt.rcParams.update({"font.family": ["Helvetica Neue", "DejaVu Sans"], "font.size": 11, "text.color": INK,
                     "axes.edgecolor": MUTED, "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": BG,
                     "xtick.labelsize": 10, "ytick.labelsize": 10, "axes.labelsize": 10.5})


def rsa_line(ax):
    ax.axvline(2048, color=GRAY, lw=1, ls=":", zorder=0)


def render(footer, name, friendly=False):
    fig = plt.figure(figsize=(7.2, 9.0), dpi=300, facecolor=BG)
    fig.text(0.06, 0.955, "How big is Shor's algorithm?", fontsize=21, fontweight="bold", va="top")
    fig.text(0.06, 0.915, "Factor 15 for real, then estimate the cost of RSA-2048", fontsize=12.5,
             color=MUTED, va="top")
    gs = fig.add_gridspec(3, 2, left=0.13, right=0.95, top=0.83, bottom=0.155, hspace=0.75, wspace=0.45)
    axes = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(6)]
    titles = (["1  Can it factor 15?", "2  What does that cost?", "3  How many logical qubits?",
               "4  How many Toffoli gates?", "5  How much protection?", "6  How many physical qubits?"]
              if friendly else
              ["1  N = 15, a = 7 (simulated)", "2  Real arithmetic vs permutation", "3  Logical qubits",
               "4  Toffoli count [GE19]", "5  Surface-code distance", "6  Physical qubits"])
    for ax, t in zip(axes, titles):
        ax.set_title(t, loc="left", fontsize=12.5, fontweight="bold", color=INK, pad=8)

    # 1 peaks of the counting register
    ax = axes[0]
    y = np.arange(2 ** A["n_count"])
    p = np.array([A["peaks"].get(int(i), 0) for i in y])
    ax.vlines(y, 0, p, color=TEAL, lw=2.4)
    ax.set_xlim(-8, 263); ax.set_ylim(0, 0.32); ax.set_xticks([0, 64, 128, 192, 255])
    ax.set_xlabel("measured y (8 qubits)"); ax.set_ylabel("probability")
    ax.text(0.98, 0.97, f"period r = {A['period']}\n15 = {A['factors'][0]} × {A['factors'][1]}", transform=ax.transAxes,
            ha="right", va="top", fontsize=10)

    # 2 gate counts: real arithmetic vs the precomputed-permutation oracle (both measured)
    ax = axes[1]
    keys, labels = ["qubits", "cx_one_multiplication", "cx"], ["qubits", "CX, one ×7", "CX, whole"]
    yy = np.arange(3)[::-1]
    for off, src, colr, lab in [(0.18, A["counts"], TEAL, "real arithmetic"), (-0.18, O["counts"], GRAY, "permutation")]:
        vals = [src[k] for k in keys]
        ax.barh(yy + off, vals, color=colr, height=0.34, label=lab)
        for v, yv in zip(vals, yy + off):
            ax.text(v * 1.25, yv, f"{v:,}", va="center", fontsize=8.5)
    ax.set_xscale("log"); ax.set_xlim(5, 1e7)
    ax.set_yticks(yy, labels); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_xlabel("count (log)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    # 3 logical qubits
    ax = axes[2]
    ax.loglog(bits, col(E3, "logical_qubits"), color=TEAL, lw=2, label="fast: 3n [GE19]")
    ax.loglog(bits, col(E3, "logical_qubits_min"), color=INK, lw=1.4, ls="--", label="fewest: 2n+3 [B03]")
    rsa_line(ax)
    ax.text(2048 * 0.92, 6189, "6,189", ha="right", va="bottom", fontsize=9.5)
    ax.set_xlabel("bits of N"); ax.set_ylabel("logical qubits")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    # 4 Toffolis
    ax = axes[3]
    ax.loglog(bits, col(E3, "toffolis"), color=TEAL, lw=2)
    rsa_line(ax)
    ax.text(2048 * 0.92, 2.6e9, "2.6 billion", ha="right", va="bottom", fontsize=9.5)
    ax.set_xlabel("bits of N"); ax.set_ylabel("Toffoli gates")

    # 5 code distance
    ax = axes[4]
    ax.step(bits, col(E3, "distance"), where="post", color=CORAL, lw=2, label="gate error 0.1 %")
    ax.step(bits, col(E4, "distance"), where="post", color=TEAL, lw=2, label="gate error 0.01 %")
    ax.set_xscale("log"); rsa_line(ax); ax.set_ylim(0, 40)
    ax.set_xlabel("bits of N"); ax.set_ylabel("code distance d")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    # 6 physical qubits + published points + today's chips
    ax = axes[5]
    ax.loglog(bits, col(E3, "physical_qubits"), color=CORAL, lw=2, label="our model, 0.1 %")
    ax.loglog(bits, col(E4, "physical_qubits"), color=TEAL, lw=2, label="our model, 0.01 %")
    ax.plot([2048], [20e6], "*", color=INK, ms=11, label="paper 2019 (20 M)")
    ax.plot([2048], [1e6], "v", color=INK, ms=7, label="preprint 2025 (< 1 M)")
    ax.axhline(1e3, color=MUTED, lw=1, ls="--")
    ax.text(9, 1.8e3, "chips today ~ 1,000", fontsize=9, color=MUTED)
    rsa_line(ax); ax.set_ylim(3e2, 3e10)
    ax.set_xlabel("bits of N"); ax.set_ylabel("physical qubits")
    ax.legend(frameon=False, fontsize=7.8, loc="upper left", handlelength=1.4)

    fig.text(0.06, 0.05, "Panels 1–2 measured in Qiskit. The permutation needs the full multiplication table: fine for 15, impossible for RSA.\n"
             "Panels 3–6 estimated: surface code, 1 µs rounds, 10 µs reaction, 1 % failure budget. RSA-2048: d = 31,\n"
             "18 M qubits, 6 h. Dotted line = RSA-2048. [B03] Beauregard 2003 · [GE19] Gidney & Ekerå 2019 · Gidney 2025 (preprint)",
             fontsize=8, color=MUTED, linespacing=1.4)
    if friendly:
        cat = load_mascot()
        w = 0.12
        cax = fig.add_axes([0.83, 0.855, w, w * (396 / 439) * (7.2 / 9.0)])
        cax.axis("off")
        if cat is not None:
            cax.imshow(cat)
        fig.text(0.75, 0.885, "Count first,\nbuild later!", fontsize=9.5, ha="center", va="center",
                 color=INK, linespacing=1.3,
                 bbox=dict(boxstyle="round,pad=0.5,rounding_size=0.9", fc="white", ec=TEAL, lw=1.5))
    fig.text(0.06, 0.02, footer, fontsize=11, color=MUTED)
    fig.text(0.95, 0.02, "@qalycat", fontsize=11, color=MUTED, ha="right")
    png = OUT / f"{name}.png"
    fig.savefig(png, facecolor=BG); plt.close(fig)
    Image.open(png).convert("RGB").resize((1080, 1350), Image.LANCZOS).save(OUT / f"{name}.jpg", quality=92)
    return png


print(render("Code: link in bio", "shor-resources-en-instagram", friendly=True))
print(render("Code: github.com/MaldoAlberto/shor-resource-estimate", "shor-resources-en-linkedin", friendly=True))
