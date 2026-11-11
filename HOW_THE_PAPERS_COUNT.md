# How the papers got "20 million qubits, 8 hours" and "under 1 million, under a week"

Two papers set the numbers everyone quotes for breaking RSA-2048 with Shor's algorithm. This page
retraces their arithmetic step by step, with the location of every number in each paper, so you can
check it and explain it.

Both papers assume the same machine: a flat grid of qubits, each talking only to its neighbours,
0.1 % error per gate, 1 µs per error-correction round, and 10 µs for the classical control to react
to a measurement. The first is peer reviewed (Quantum, 2021); the second is a 2025 preprint.

---

## Paper 1: Gidney & Ekerå 2019, "20 million qubits, 8 hours"

C. Gidney, M. Ekerå, Quantum 5, 433 (2021), arXiv:1905.09749. Sections 2.10–2.14, Tables 2 and 3.

### Step 1. Shrink the algorithm before counting anything

The cost of Shor's algorithm is almost all modular exponentiation: a long chain of multiplications.
The paper stacks known tricks to make that chain cheaper:

| Trick | Plain words |
|---|---|
| Ekerå–Håstad variant | Turns factoring into a "short discrete logarithm", so the exponent has about 1.5n bits (n_e = 3029 for n = 2048) instead of 2n. |
| Windowed arithmetic | Looks up precomputed values in a table, so 5 × 5 bits are handled in one step instead of one bit at a time. |
| Coset representation | Uses an approximate form of modular arithmetic that skips the expensive "reduce mod N" step, at a tiny controlled error. |
| Oblivious carry runways | Cuts the 2048-bit addition into pieces of 1024 bits that run in parallel. |

The result is one repeated building block, the lookup addition. RSA-2048 needs about
0.1 · n_e · n = 0.1 · 3029 · 2048 ≈ **620,000 lookup additions** (Eq. 1). The abstract gives the
totals: 3n + 0.002 n lg n logical qubits (6,189) and 0.3n³ + 0.0005 n³ lg n Toffolis (2.6 billion).

### Step 2. Draw the chip layout, then count patches

Each logical qubit is a square surface-code patch. The paper draws the whole machine as a grid of
patches (Section 2.10, Figures 5 and 6):

- Each of the 2 parallel pieces is 113 patches wide. 14 CCZ factories of 15 × 8 patches sit in 2 rows of 7, with gaps for routing.
- Each piece is 63 rows tall: 16 rows of factories, 30 rows of register data, 8 for routing, 6 for CCZ fix-ups and 3 for the adder.
- Total: 113 × 2 × 63 = **14,238 patches**.

### Step 3. Patches to physical qubits

"A logical qubit covers 2(d+1)² physical qubits" (Section 2.14). At code distance d = 27 that is
2 · 28² = **1,568 physical qubits per patch**, so 14,238 × 1,568 ≈ **22.3 million**. The paper says
"approximately 23 million" and explains that its script tunes the layout further (the level-1
factory distance drops from 17 to 15) and lands on **20 million** (Table 3).

Where the 22.3 million go, by rows of the layout:

| Region | Rows | Physical qubits |
|---|---|---|
| Register data | 30 | 10.6 M |
| Magic-state factories | 16 | 5.7 M |
| Routing | 8 | 2.8 M |
| CCZ fix-ups | 6 | 2.1 M |
| Adder | 3 | 1.1 M |

### Step 4. Runtime: one step times the number of steps

One lookup addition takes about **37 ms** (Section 2.11):

- Lookup: limited by code distance: 1 µs × d/2 × 2^(c_exp + c_mul) ≈ 14 ms.
- Addition: limited by the 10 µs reaction time, one carry at a time: 2(c_sep + c_pad) × 10 µs ≈ 22 ms.
- Everything else: about 1 ms.

620,000 × 37 ms ≈ **6.4 hours** (Eq. 5, "approximately 7 hours"). The tuned script gives
**5.1 hours per run** (Table 3).

### Step 5. Accept some risk, then pay for retries

The paper does not demand a perfect run. At n = 2048 it accepts a **31 % chance that a run fails**
(Table 3, "retry risk"). A failed run is detected and repeated, so the expected time is
5.1 h / (1 − 0.31) ≈ **7.4 hours**. The title rounds this up to 8 hours. A single successful run
is enough to factor with at least 99 % success in classical post-processing (Section 3).

This choice explains why the paper uses d = 27, while our simpler model picks d = 31: our model asks
for a 1 % failure budget, not 31 %.

---

## Paper 2: Gidney 2025, "under 1 million qubits, under a week"

C. Gidney, arXiv:2505.15917 (2025, preprint, not peer reviewed). Section 3.2 and Table 5.

### Step 1. Use far fewer logical qubits, pay with more gates

- **Approximate residue arithmetic** (Chevignard, Fouque, Schrottenloher 2024) does the arithmetic modulo many small primes instead of on the full 2048-bit number. This needs far fewer logical qubits: about 1,400 at the peak (1,409 in Section 3.2; Table 5 lists 1,399), down from 6,189.
- The price is more work: **6.5 billion Toffolis** (Table 5), more than twice the 2019 design.

### Step 2. Store idle qubits more densely

Most logical qubits sit idle most of the time. **Yoked surface codes** store idle ("cold") qubits at
**430 physical qubits** each, instead of **1,352** for an active ("hot") patch at d = 25, using
2(d+1)² = 2 · 26² (Section 3.2). That is roughly 3× denser.

### Step 3. Smaller magic-state factories

**Magic-state cultivation** (Gidney, Shutty, Jones 2024) replaces the first distillation stage. The
factories shrink, and the design uses only **6 factories** instead of 28.

### Step 4. Add up the regions

| Region | Count | Physical qubits each | Total |
|---|---|---|---|
| Cold storage | 1,280 logical qubits | 430 | 550,400 |
| Hot storage | 131 logical qubits | 1,352 | 177,112 |
| Compute (6 factories + workspace) | 7 × 18 = 126 patches | 1,352 | 170,352 |
| **Total** | | | **897,864** |

The paper rounds up "for slack" and reports **1 million**.

### Step 5. Runtime: per shot, times shots, divided by success

- One shot takes **12.07 hours**. Six factories give one CCZ state every 25 µs, which sets the pace.
- The algorithm needs several shots: **9.2 on average** (Table 5).
- Each shot finishes without a logical error with probability **93.3 %**. The target is 10⁻¹⁵ error per logical qubit per round, over about 6.9 × 10¹³ qubit-rounds.
- 12.07 h × 9.2 / 0.933 = **4.96 days**, reported as "less than a week".

---

## The two papers side by side

| | 2019 | 2025 |
|---|---|---|
| Logical qubits | 6,189 (formula) | about 1,400 at peak |
| Toffolis | 2.6 billion | 6.5 billion |
| Factories | 28 | 6 |
| Physical qubits per idle logical qubit | 1,568 (d = 27) | 430 (yoked, cold) |
| Physical qubits | 20 million | 0.9 million, reported as 1 million |
| Expected time | 7.4 h, title "8 hours" | 4.96 days, "less than a week" |

The 20× drop in qubits came from algorithms and error-correction design, not from a better chip:
both papers assume the same 0.1 % gate error. The trade-off is time: 2025 does more Toffolis with
fewer factories, so it runs about 16× longer.

## What today's hardware shows

Λ is the factor by which the logical error drops each time d grows by 2. With the threshold
heuristic Λ ≈ p_th / p, the papers' 0.1 % gate error and a 1 % threshold correspond to Λ ≈ 10.

| Group, platform | Code | Λ | Source |
|---|---|---|---|
| USTC, superconducting | surface, d = 7 | 1.40 ± 0.06 | PRL 135, 260601 (2025) |
| Google, superconducting, 105 qubits | surface, d = 7 | 2.14 ± 0.02 | Nature 638, 920 (2025) |
| Harvard/MIT/QuEra, 448 neutral atoms | surface, d = 3 → 5 | 2.14 ± 0.13 | Nature (2025), arXiv:2506.20661 |
| Quantinuum, trapped ions | color, d = 7 | 12 (+3/−2) | arXiv:2610.24077 (preprint) |

These values come from different codes, protocols and decoders, so they are not a league table.
The ion result reaches a similar suppression on a different code (color code), decoded offline, on
73 qubits, not yet at speed or scale. The superconducting surface-code chips have the speed but
Λ ≈ 1.4–2.1. The papers need surface-code Λ ≈ 10, 1 µs rounds and about a million qubits together.
That combination is the hardware gap.

## What came after 2025

2026 preprints (not peer reviewed) lower the qubit count further, each with its own machine
assumptions:
- Webster et al. (Pinnacle, arXiv:2602.11457): fewer than 100,000 qubits with qLDPC codes at 0.1 %
  error, 1 µs cycle and 10 µs reaction, but qLDPC codes need longer-range links than a flat
  nearest-neighbour grid.
- Mundada et al. (arXiv:2604.06319): 381,000 qubits and 9.2 days with a heterogeneous design; the
  abstract does not state the error rate or cycle time.
- Xue & Covey (arXiv:2605.03951): half a million qubits on a modular atomic processor, assuming
  1 ms measurements instead of 1 µs rounds.

## How to explain it in one minute

1. **Shrink the recipe.** Mathematical tricks turn the algorithm into about 620,000 copies of one step.
2. **Draw the factory floor.** Place every logical qubit and every magic-state factory on a grid, then count the squares.
3. **Price each square.** One reliable square costs about 1,500 noisy qubits in the 2019 design, because the noisy qubits check each other.
4. **Time one step, multiply.** 37 ms × 620,000 ≈ 6 hours. Then add retries for runs that fail.
5. **Ask what changed.** The 2025 paper uses fewer squares and cheaper storage for idle ones, and it accepts a longer run.
