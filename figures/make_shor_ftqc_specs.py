"""Figure "What quantum computer does Shor need?": the FTQC spec sheet behind the RSA-2048 estimate
and how the answer moves when each spec moves. Companion to make_shor_resources.py.
Data: ../shor_estimate.py (model + published formulas; nothing here is hardware data
except measured Lambda values: surface code USTC 1.40, Google 2.14, Harvard 2.14; color code Quantinuum 12 (preprint)).
Cream background = intermediate category. English, LinkedIn + Instagram (Qaly, "link in bio").
Usage (from the repo root): python3 figures/make_shor_ftqc_specs.py <output folder>
"""
import pathlib, sys

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

N_BITS = 2048
E = se.estimate(N_BITS, 1e-3)
D_RSA = E["distance"]
ROUNDS = E["meas_depth"] * se.REACTION_ROUNDS                                # reaction-limited: 10 rounds per layer
NEED = se.FAIL_BUDGET / (E["logical_qubits"] * ROUNDS)                     # p_L budget per qubit per round

plt.rcParams.update({"font.family": ["Helvetica Neue", "DejaVu Sans"], "font.size": 11, "text.color": INK,
                     "axes.edgecolor": MUTED, "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": BG,
                     "xtick.labelsize": 10, "ytick.labelsize": 10, "axes.labelsize": 10.5})


def render(footer, name, friendly=True):
    fig = plt.figure(figsize=(7.2, 9.0), dpi=300, facecolor=BG)
    fig.text(0.06, 0.955, "What quantum computer does Shor need?", fontsize=20, fontweight="bold", va="top")
    fig.text(0.06, 0.915, "Our model's specs for RSA-2048, and if they slip", fontsize=12,
             color=MUTED, va="top")
    gs = fig.add_gridspec(3, 2, left=0.13, right=0.95, top=0.80, bottom=0.155, hspace=0.75, wspace=0.45)
    axes = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(6)]
    titles = ["1  Our model's spec sheet", "2  Why stay below threshold?", "3  Price of one logical qubit",
              "4  What if gates are noisier?", "5  What if the clock is slower?", "6  Where do the qubits go?"]
    for ax, t in zip(axes, titles):
        ax.set_title(t, loc="left", fontsize=12.5, fontweight="bold", color=INK, pad=8)

    # 1 spec sheet
    ax = axes[0]; ax.axis("off")
    specs = [("gate error p", "0.1 %"), ("threshold", "1 %"), ("code", "surface, 2D grid"),
             ("round time", "1 µs"), ("decoder reacts in", "10 µs"), ("run may fail (ours)", "1 %"),
             ("routing (ours)", "× 1.5"), ("factories (ours)", "4 × 12")]
    for i, (k, v) in enumerate(specs):
        yy = 1.0 - i * 0.135
        ax.text(0.0, yy, k, fontsize=10, color=MUTED, transform=ax.transAxes, va="top")
        ax.text(1.0, yy, v, fontsize=10, color=INK, fontweight="bold", transform=ax.transAxes, va="top", ha="right")

    # 2 logical error vs distance
    ax = axes[1]
    d = np.arange(3, 52, 2)
    for p, c, lab, x, dy, ha in [(1e-4, TEAL, "p = 0.01 %", 4, 1e-17, "left"), (1e-3, INK, "p = 0.1 %", 42, 1e-3, "right"),
                                 (5e-3, CORAL, "p = 0.5 %\n(Λ = 2 in model)", 51, 1e4, "right")]:
        ax.semilogy(d, se.p_logical(d, p), color=c, lw=2)
        ax.text(x, se.p_logical(x, p) * dy, lab, color=c, fontsize=8.5, ha=ha, va="top" if dy < 1 else "bottom")
    ax.axhline(NEED, color=MUTED, lw=1, ls="--")
    ax.text(51, NEED * 3, "RSA-2048 needs", fontsize=8.5, color=MUTED, ha="right", va="bottom")
    ax.plot([D_RSA], [se.p_logical(D_RSA, 1e-3)], "o", color=INK, ms=6)
    ax.set_ylim(1e-30, 1); ax.set_yticks([1, 1e-10, 1e-20, 1e-30])
    ax.set_xlabel("code distance d"); ax.set_ylabel("logical error / round")

    # 3 physical qubits per logical qubit
    ax = axes[2]
    ax.plot(d, 2 * d ** 2, color=TEAL, lw=2)
    ax.plot([D_RSA], [2 * D_RSA ** 2], "o", color=INK, ms=6)
    ax.plot([27], [1568], "s", color=CORAL, ms=6)
    ax.text(4, 5800, f"our model, d = {D_RSA}: 2·{D_RSA}² = {2 * D_RSA ** 2:,}\n(×1.5 routing ≈ {2 * D_RSA ** 2 * se.ROUTING:,.0f})\n"
            "paper 2019, d = 27: 2·28² = 1,568", fontsize=8.5, va="top", linespacing=1.3)
    ax.set_xlabel("code distance d"); ax.set_ylabel("physical per logical")
    ax.set_ylim(0, 6000)

    # 4 physical qubits vs gate error
    ax = axes[3]
    ps = np.logspace(-4, np.log10(7.5e-3), 50)
    q = [se.estimate(N_BITS, p) for p in ps]
    ok = [(p, e["physical_qubits"]) for p, e in zip(ps, q) if e]
    ax.loglog(*zip(*ok), color=TEAL, lw=2)
    ax.axvspan(ok[-1][0], 1e-2, color=CORAL, alpha=0.15, lw=0)
    ax.text(0.99, 0.04, "no fit\nd ≤ 201", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=CORAL)
    ax.plot([1e-3], [E["physical_qubits"]], "o", color=INK, ms=6)
    ax.text(1e-3, E["physical_qubits"] * 2.2, "18 M", ha="center", fontsize=9.5)
    E5 = se.estimate(N_BITS, 5e-3)
    ax.plot([5e-3], [E5["physical_qubits"]], "o", color=CORAL, ms=6)
    ax.text(4.3e-3, E5["physical_qubits"] * 1.6, f"{E5['physical_qubits'] / 1e6:.0f} M, d = {E5['distance']}", ha="right", va="bottom", fontsize=8.5, color=CORAL)
    ax.set_xlim(1e-4, 1e-2); ax.set_ylim(1e6, 1e9)
    ax.set_xlabel("gate error p"); ax.set_ylabel("physical qubits")

    # 5 runtime vs cycle time
    ax = axes[4]
    cyc = np.logspace(-1, 3, 50)                      # microseconds
    days = ROUNDS * cyc / 3.6e9 / 24
    ax.loglog(cyc, days, color=TEAL, lw=2)
    for c, lab, ha in [(1, "1 µs\n6 hours", "left"), (1000, "1 ms\n~8 months", "right")]:
        v = ROUNDS * c / 3.6e9 / 24
        ax.plot([c], [v], "o", color=INK, ms=6)
        ax.text(c * (1.6 if ha == "left" else 0.8), v * (1 if ha == "left" else 1.6), lab, fontsize=9, ha=ha,
                va="top" if ha == "left" else "bottom")
    ax.set_ylim(top=3e3); ax.set_xlabel("round time (µs), reaction = 10 rounds"); ax.set_ylabel("days")

    # 6 where the qubits go
    ax = axes[5]
    data = E["logical_qubits"] * 2 * D_RSA ** 2
    routing = data * (se.ROUTING - 1)
    fact = E["physical_qubits"] - data - routing
    labels, vals = ["data", "routing", "factories*"], [data, routing, fact]
    yy = np.arange(3)[::-1]
    ax.barh(yy, vals, color=[TEAL, GRAY, CORAL], height=0.6)
    ax.set_yticks(yy, labels); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 1.6e7); ax.set_xticks([0, 5e6, 1e7, 1.5e7], ["0", "5 M", "10 M", "15 M"])
    for v, yv in zip(vals, yy):
        ax.text(v + 3e5, yv, f"{v / 1e6:.1f} M" if v > 1e6 else f"{v / 1e3:.0f} k", va="center", fontsize=10)
    ax.set_xlabel("physical qubits (p = 0.1 %)")

    fig.text(0.06, 0.05, "Model: p_L = 0.1 (p / 1 %)^((d+1)/2) (Gidney & Ekerå 2019), 2d² qubits per patch, ×1.5 routing, Λ = 10.\n"
             "Measured surface-code Λ: 1.4 USTC, 2.14 Google, 2.14 Harvard. Color code: 12 Quantinuum (preprint,\n"
             "73 qubits, offline decoding). *Factories: crude guess (about 25 % of the chip in the 2019 layout).",
             fontsize=8, color=MUTED, linespacing=1.4)
    if friendly:
        cat = load_mascot()
        w = 0.10
        cax = fig.add_axes([0.85, 0.85, w, w * (396 / 439) * (7.2 / 9.0)])
        cax.axis("off")
        if cat is not None:
            cax.imshow(cat)
        fig.text(0.86, 0.90, "Specs first,\nthen the count!", fontsize=9, ha="right", va="center",
                 color=INK, linespacing=1.3,
                 bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.9", fc="white", ec=TEAL, lw=1.5))
    fig.text(0.06, 0.02, footer, fontsize=11, color=MUTED)
    fig.text(0.95, 0.02, "@qalycat", fontsize=11, color=MUTED, ha="right")
    png = OUT / f"{name}.png"
    fig.savefig(png, facecolor=BG); plt.close(fig)
    Image.open(png).convert("RGB").resize((1080, 1350), Image.LANCZOS).save(OUT / f"{name}.jpg", quality=92)
    return png


print(render("Code: link in bio", "shor-ftqc-specs-en-instagram"))
print(render("Code: github.com/MaldoAlberto/shor-resource-estimate", "shor-ftqc-specs-en-linkedin"))
