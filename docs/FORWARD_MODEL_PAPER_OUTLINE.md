# Forward-model paper outline

**Audience:** first draft structure for the October 2026 manuscript  
**Paper type:** methods / software  
**Claim focus:** integrated, validated morphological forward platform—not novelty of dilation itself  
**Out of scope for this paper:** inverse reconstruction (future work)

One-sentence claim to keep the draft on track:

> We present a transparent, reproducible morphological AFM forward simulator and show that it matches analytical tip geometry, SEM-derived tip profiles, and real AFM measurements of a known sample within quantified residuals.

---

## Recommended structure

### Title

Something like: *An interactive morphological forward model for AFM tip–sample convolution, with geometry and experimental validation*

### Abstract (150–250 words)

Write this last. Cover, in order: problem → what you built → how you validated → main quantitative results → why it matters. One sentence each is enough.

### 1. Introduction

- AFM images mix sample and tip geometry (broadening, inaccessible valleys, double-tip ghosts).
- Why geometric forward models still matter for teaching, debugging, and later inverse work.
- Gap: most tools are either opaque, hard to share, or not validated end-to-end against real data.
- Contribution list (3–4 bullets):
  - browser + Python forward model;
  - ideal / non-ideal / OBJ probes;
  - mesh / SEM / `.spm` validation;
  - export / reproducibility path toward inverse benchmarks.
- One sentence on scope: contact-mode height-field morphology only; inverse is future work.

### 2. Related work

Keep short (1–1.5 pages). Group by:

- classical tip convolution / morphology (Villarrubia);
- AFM simulation / analysis software (e.g. Gwyddion);
- tip estimation and reconstruction (classical + recent differentiable/ML work—only enough to position this paper as forward, not inverse);
- interactive teaching / visualization tools if relevant.

End with: *what prior work does not provide that yours does*.

### 3. Methods

This is the core. Subsections:

#### 3.1 Geometric forward model

Dilation equation, erosion as optional dual, height-field assumptions, what you deliberately exclude (deformation, feedback, overhangs).

#### 3.2 Probe and surface representations

Analytical tips; faceted/non-ideal; OBJ → lower/upper envelope; synthetic surfaces used in the paper.

#### 3.3 Implementation

2D vs 3D simulators; Python mirror; grid (80×80 / 200 nm); exports (CSV, PNG, OBJ, JSON). One diagram of the pipeline is enough.

#### 3.4 Validation protocol

Three tiers, clearly separated:

1. analytical ↔ exported mesh;
2. SEM probe profile ↔ simulator probe;
3. forward prediction ↔ real `.spm` on a known sample.

State metrics (RMSE, max abs error, coverage), tolerances, and how scale/alignment were handled.

### 4. Results

Lead with evidence, not features.

#### 4.1 Geometry consistency

Cone / pyramid mesh validation numbers already available.

#### 4.2 Probe-profile comparison

SEM extraction vs simulator probe (residuals + uncertainty).

#### 4.3 Forward-model vs real AFM

Predicted map vs `.spm`; residual maps; line profiles.

#### 4.4 Cross-implementation checks

Browser ↔ Python agreement (if included in the frozen package).

#### 4.5 Negative controls / resolution

Only if included in the frozen package.

Each subsection: 1–2 figures + a short table of numbers.

### 5. Discussion

- What the validation supports vs what it does not.
- Sources of residual error (discretization, SEM extraction, tip kernel footprint, real AFM dynamics not modeled).
- Use cases: teaching, probe artifact intuition, synthetic data for later inverse work.
- Limitations paragraph (honest and specific).
- One short paragraph pointing to the inverse roadmap—without overselling.

### 6. Conclusions

3–5 sentences: what was built, what was validated, main quantitative takeaway, next step (inverse on known samples).

### Availability

Software version, license, how to run validation commands, where data/figures live.

### References

Existing short list plus whatever the venue requires.

---

## Figures to plan early (write around these)

| Fig | Content |
|-----|---------|
| 1 | Pipeline: surface + tip → dilation → measured map |
| 2 | Interactive UI screenshot (optional; one is enough) |
| 3 | Mesh vs analytical profiles (cone / pyramid) |
| 4 | SEM tip profile vs simulator tip |
| 5 | True sample / simulated AFM / real `.spm` / residual |
| 6 | Line-profile overlays for the real-data case |

Plus one results table of RMSE / max error / coverage.

---

## Writing order for a first paper

Do **not** start with the Introduction. Order that usually works:

1. Methods 3.1–3.3 (already known)
2. Results figures + captions
3. Methods 3.4 (validation protocol matching the figures)
4. Discussion + Limitations
5. Introduction + Related work
6. Abstract + Title
