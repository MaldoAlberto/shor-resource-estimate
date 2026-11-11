"""Beginner figure v2 "Can a quantum computer break RSA yet?" as an explainer: each row has a plain
explanation (what we did), a chart (the evidence) and a takeaway (what it means). Jargon is defined
in the row where it first appears. Same verified data as make_shor_beginner.py.
Blue = beginner category. Usage (from the repo root): python3 figures/make_shor_beginner_v2.py <output folder>
"""
import json, pathlib, sys, textwrap

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
A = json.loads((REPO / "arith.json").read_text())
A["peaks"] = {int(k): v for k, v in A["peaks"].items()}

INK, MUTED = "#13203A", "#5B6475"
TEAL, CORAL, GRAY = "#0A9396", "#C94A45", "#B7BDC8"
BG = "#EAF1FA"  # blue = beginner category

plt.rcParams.update({"font.family": ["Helvetica Neue", "DejaVu Sans"], "font.size": 10, "text.color": INK,
                     "axes.edgecolor": MUTED, "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": BG,
                     "xtick.labelsize": 9, "ytick.labelsize": 9.5, "axes.labelsize": 9.5})

ROWS = [
    {"q": "1  Does it work on a small number?",
     "body": "Shor's algorithm finds the factors of a number by finding a pattern that repeats. "
             "We built the full circuit for 15 from real adders and simulated it on a laptop.",
     "take": "Four evenly spaced peaks: the pattern repeats every 4 steps, and that gives 15 = 3 × 5."},
    {"q": "2  How many qubits for RSA-2048?",
     "body": "RSA-2048 protects much of the web. Qubits make mistakes, so one reliable qubit "
             "is a team of about 1,500 noisy ones checking each other (2019 design). "
             "The 2025 design packs idle ones into about 430 each.",
     "take": "Better designs cut the need 20× since 2019, but chips are still 100× too small or more."},
    {"q": "3  How long would it run?",
     "body": "Even with enough qubits, the job is billions of steps, each waiting for "
             "errors to be found and fixed.",
     "take": "Hours to a week. Time is not the barrier; the size of the machine is."},
    {"q": "4  How fast do errors drop?",
     "body": "Error correction adds layers of protection. For these estimates, each extra "
             "layer must cut errors about 10×.",
     "take": "Fast chips reach about 2×. Ions reach 12× on a different code, on 73 qubits, "
             "not yet at speed or scale. No machine has all three together yet."},
]


def hbars(ax, labels, vals, colors, fmt, log=True):
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=colors, height=0.62)
    ax.set_yticks(y, labels); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    if log:
        ax.set_xscale("log")
    for v, yy in zip(vals, y):
        ax.text(v * 1.4 if log else v + 0.3, yy, fmt(v), va="center", fontsize=9.5)


def chart(i, ax):
    if i == 0:
        y = np.arange(2 ** A["n_count"])
        ax.vlines(y, 0, [A["peaks"].get(int(k), 0) for k in y], color=TEAL, lw=2.4)
        ax.set_xlim(-8, 263); ax.set_ylim(0, 0.3); ax.set_xticks([0, 64, 128, 192, 255]); ax.set_yticks([0, 0.25])
        ax.set_ylabel("chance"); ax.set_xlabel("measured result")
    elif i == 1:
        hbars(ax, ["paper 2019", "preprint 2025", "preprint 2026", "largest chips"], [20e6, 1e6, 1e5, 1e3],
              [INK, INK, INK, GRAY], lambda v: {20e6: "20 M", 1e6: "< 1 M", 1e5: "< 100 k", 1e3: "~1,000"}[v])
        ax.set_xlim(3e2, 3e9); ax.set_xticks([1e3, 1e6, 1e9], ["thousand", "million", "billion"]); ax.minorticks_off()
    elif i == 2:
        hbars(ax, ["paper 2019", "preprint 2025"], [8, 168], [INK, INK],
              lambda v: {8: "8 hours", 168: "< 1 week"}[v])
        ax.set_xlim(1, 3e3); ax.set_xticks([1, 24, 168], ["1 h", "1 day", "1 week"]); ax.minorticks_off()
    else:
        labs = ["USTC chip", "Google chip", "Harvard atoms", "ions, color code*", "needed"]
        vals = [1.40, 2.14, 2.14, 12, 10]
        hbars(ax, labs, vals, [TEAL] * 4 + [CORAL], lambda v: f"÷{v:g}", log=False)
        ax.set_xlim(0, 15); ax.set_xticks([0, 5, 10])
        ax.set_xlabel("errors cut per extra layer")


def render(footer, name):
    fig = plt.figure(figsize=(7.2, 9.0), dpi=300, facecolor=BG)
    fig.text(0.06, 0.96, "Can a quantum computer\nbreak RSA yet?", fontsize=21, fontweight="bold", va="top",
             linespacing=1.15)
    fig.text(0.06, 0.875, "Four questions, answered with data", fontsize=12, color=MUTED, va="top")

    tops = [0.83, 0.66, 0.46, 0.30]           # row 2 has the longest text
    for i, r in enumerate(ROWS):
        y0 = tops[i]
        fig.text(0.06, y0, r["q"], fontsize=12.5, fontweight="bold", va="top")
        body = textwrap.fill(r["body"], 46)
        fig.text(0.06, y0 - 0.03, body, fontsize=9.4, color=INK, va="top", linespacing=1.35)
        yt = y0 - 0.03 - 0.0195 * (body.count("\n") + 1) - 0.008
        fig.text(0.06, yt, textwrap.fill(r["take"], 48), fontsize=9.2, color=TEAL, va="top",
                 linespacing=1.3, fontweight="bold")
        ax = fig.add_axes([0.66 if i == 3 else (0.62 if i else 0.55), y0 - 0.145, 0.29 if i == 3 else (0.33 if i else 0.40), 0.11])
        chart(i, ax)

    fig.text(0.06, 0.098, "Not yet. And now we can say exactly what is missing.", fontsize=12, fontweight="bold")
    fig.text(0.06, 0.037,
             "Chips and atoms: surface code, as the estimates assume. Preprint 2026: qLDPC code, needs long-range links.\n"
             "*Quantinuum ions: color code, preprint, 73 qubits, decoded offline, rounds far slower than 1 µs.\n"
             "Gidney & Ekerå 2019 · Gidney 2025 (preprint) · Webster et al. 2026 (preprint) · He et al. PRL 2025 · Google Nature 2025 ·\n"
             "Bluvstein et al. Nature 2025 · Quantinuum arXiv:2610.24077. Circuit for 15: our code, 18 qubits.",
             fontsize=7, color=MUTED, linespacing=1.4)
    cat = load_mascot()
    w = 0.12
    cax = fig.add_axes([0.83, 0.865, w, w * (396 / 439) * (7.2 / 9.0)])
    cax.axis("off")
    if cat is not None:
        cax.imshow(cat)
    fig.text(0.81, 0.925, "Let's count\nthe gap!", fontsize=9.5, ha="right", va="center", color=INK, linespacing=1.3,
             bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.9", fc="white", ec=TEAL, lw=1.5))
    fig.text(0.06, 0.015, footer, fontsize=10.5, color=MUTED)
    fig.text(0.95, 0.015, "@qalycat", fontsize=10.5, color=MUTED, ha="right")
    png = OUT / f"{name}.png"
    fig.savefig(png, facecolor=BG); plt.close(fig)
    Image.open(png).convert("RGB").resize((1080, 1350), Image.LANCZOS).save(OUT / f"{name}.jpg", quality=92)
    return png


print(render("Code: link in bio", "shor-beginner-v2-en-instagram"))
print(render("Code: github.com/MaldoAlberto/shor-resource-estimate", "shor-beginner-v2-en-linkedin"))
