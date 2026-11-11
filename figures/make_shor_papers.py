"""Beginner figure "How did they get 20 million and 8 hours?": retraces the arithmetic of
Gidney & Ekerå 2019 (Quantum 5, 433; Sec. 2.10-2.14, Tables 2-3) and Gidney 2025 (arXiv:2505.15917,
Sec. 3.2, Table 5). Every number is quoted from, or multiplied out of, those papers.
Explainer: ../HOW_THE_PAPERS_COUNT.md. Blue = beginner category.
Usage (from the repo root): python3 figures/make_shor_papers.py <output folder>
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


OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "out").resolve()
OUT.mkdir(parents=True, exist_ok=True)

INK, MUTED = "#13203A", "#5B6475"
TEAL, CORAL, GRAY, SAND = "#0A9396", "#C94A45", "#B7BDC8", "#E0A458"
BG = "#EAF1FA"  # blue = beginner category

# ---------------------------------------------------------------- data from the papers
PER_ROW = 113 * 2 * 1568                       # GE19: 113 patches wide x 2 pieces x 2(27+1)^2
ROWS19 = [("stored numbers", 30), ("factories", 16), ("routing", 8), ("other", 6), ("adder", 3)]
QB25 = [("cold storage", 1280 * 430), ("hot storage", 131 * 1352), ("factories +\nwork space", 126 * 1352)]
T19 = [("620k steps × 37 ms", 620_339 * 0.037 / 3600), ("tuned, 1 run", 5.1), ("+ retries", 5.1 / 0.69)]
T25 = [("1 try", 12.07), ("× 9.2 tries", 12.07 * 9.2), ("+ failed runs", 12.07 * 9.2 / 0.933)]

plt.rcParams.update({"font.family": ["Helvetica Neue", "DejaVu Sans"], "font.size": 12, "text.color": INK,
                     "axes.edgecolor": MUTED, "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": BG,
                     "xtick.labelsize": 10.5, "ytick.labelsize": 10.5, "axes.labelsize": 11})


def hbars(ax, labels, vals, colors, fmt, xmax):
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=colors, height=0.62)
    ax.set_yticks(y, labels); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, xmax)
    for v, yy in zip(vals, y):
        ax.text(v + xmax * 0.02, yy, fmt(v), va="center", fontsize=10.5)


def render(footer, name):
    fig = plt.figure(figsize=(7.2, 9.0), dpi=300, facecolor=BG)
    fig.text(0.06, 0.955, "How did the papers\nget their numbers?", fontsize=22, fontweight="bold", va="top",
             linespacing=1.15)
    fig.text(0.06, 0.865, "Retracing the arithmetic of the two RSA-2048 estimates", fontsize=12.5, color=MUTED, va="top")
    gs = fig.add_gridspec(2, 2, left=0.23, right=0.96, top=0.78, bottom=0.23, hspace=0.6, wspace=0.95)
    axes = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]
    titles = ["2019 paper: where\nthe qubits go", "2025 preprint: where\nthe qubits go", "2019: how long?", "2025: how long?"]
    for ax, t in zip(axes, titles):
        ax.set_title(t, loc="left", fontsize=13, fontweight="bold", color=INK, pad=10, x=-0.55)

    hbars(axes[0], [k for k, _ in ROWS19], [r * PER_ROW / 1e6 for _, r in ROWS19],
          [TEAL, CORAL, GRAY, GRAY, GRAY], lambda v: f"{v:.1f} M", 14)
    axes[0].set_xlabel("million physical qubits")
    axes[0].text(0.98, 0.02, "sum 22 M\nscript: 20 M", transform=axes[0].transAxes, ha="right", va="bottom",
                 fontsize=10, color=INK)

    hbars(axes[1], [k for k, _ in QB25], [v / 1e3 for _, v in QB25], [TEAL, SAND, CORAL],
          lambda v: f"{v:.0f} k", 800)
    axes[1].set_xlabel("thousand physical qubits")
    axes[1].text(0.98, 0.02, "sum 898 k\nreported: 1 M", transform=axes[1].transAxes, ha="right", va="bottom",
                 fontsize=10, color=INK)

    hbars(axes[2], [k for k, _ in T19], [v for _, v in T19], [GRAY, GRAY, TEAL], lambda v: f"{v:.1f} h", 10)
    axes[2].set_xlabel("hours  (title: \"8 hours\")"); axes[2].set_xticks([0, 4, 8])

    hbars(axes[3], [k for k, _ in T25], [v / 24 for _, v in T25], [GRAY, GRAY, TEAL], lambda v: f"{v:.1f} d", 7)
    axes[3].set_xlabel("days  (\"less than a week\")")

    fig.text(0.06, 0.075,
             "Recipe: draw the chip, count the squares, price each square (2019: 1,568 noisy qubits),\n"
             "time one step, add retries. 2025: fewer squares, idle ones stored 3× denser, smaller factories.\n"
             "Same machine in both: 0.1 % gate error, 1 µs rounds, 10 µs reaction time. 2019 run fails 31 % → retry.\n"
             "Gidney & Ekerå, Quantum 2021 (Sec. 2.10–2.14, Table 3) · Gidney, arXiv:2505.15917 (Sec. 3.2)",
             fontsize=8, color=MUTED, linespacing=1.45)
    cat = load_mascot()
    w = 0.13
    cax = fig.add_axes([0.82, 0.855, w, w * (396 / 439) * (7.2 / 9.0)])
    cax.axis("off")
    if cat is not None:
        cax.imshow(cat)
    fig.text(0.80, 0.915, "Draw the chip,\ncount the squares,\ntime one step!", fontsize=9.5, ha="right",
             va="center", color=INK, linespacing=1.3,
             bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.9", fc="white", ec=TEAL, lw=1.5))
    fig.text(0.06, 0.02, footer, fontsize=11, color=MUTED)
    fig.text(0.95, 0.02, "@qalycat", fontsize=11, color=MUTED, ha="right")
    png = OUT / f"{name}.png"
    fig.savefig(png, facecolor=BG); plt.close(fig)
    Image.open(png).convert("RGB").resize((1080, 1350), Image.LANCZOS).save(OUT / f"{name}.jpg", quality=92)
    return png


print(render("Code: link in bio", "shor-papers-en-instagram"))
print(render("Code: github.com/MaldoAlberto/shor-resource-estimate", "shor-papers-en-linkedin"))
