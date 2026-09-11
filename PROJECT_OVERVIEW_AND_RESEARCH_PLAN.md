# AFM Forward Simulation, Validation, and Inverse Reconstruction Roadmap

**Project status:** Active research and software development  
**Current emphasis:** Forward-model manuscript draft by 1 October 2026; submission readiness by 31 October 2026; inverse-problem foundation by 27 November 2026  
**Capacity:** Approximately 10 hours per week (about 8 focused research hours plus 2 hours for meetings, review, documentation, and contingency)  


## Executive summary

Atomic force microscopy (AFM) measures surface topography by moving a physical probe across a sample. Because the probe has a finite shape, the recorded image is not a direct copy of the surface: it combines sample geometry with probe geometry. Narrow gaps may be inaccessible, sharp features may appear broadened or rounded, and damaged or multi-apex probes can create duplicated features.

This collaborative project develops a transparent, browser-based platform for exploring and quantifying these effects. Its core is a contact-mode height-field forward model based on mathematical morphology. Given a known surface and probe geometry, the software computes the AFM image that the probe would measure. The platform currently provides:

- interactive 2D line-scan and 3D raster-scan simulators;
- parameterized ideal, non-ideal, asymmetric, multi-apex, concave, and faceted probe models;
- synthetic surfaces spanning simple calibration shapes and more complex topographies;
- import and export of Wavefront OBJ geometry;
- export of true and measured height maps, metadata, and error statistics;
- a Python implementation for repeatable simulation and future batch generation; and
- independent geometry-validation tooling that compares exported meshes with the analytical equations used by the simulator.

The immediate research objective is to establish a reliable and reproducible forward model and complete a forward-model manuscript this semester (complete draft by 1 October 2026; submission-ready by 31 October 2026). Forward-model validation uses real AFM measurement data (Bruker `.spm` files) from known calibration samples: the simulator predicts what the measurement should look like given an independently characterized probe, and the prediction is compared against the real measurement. The longer-term objective is the inverse problem: recovering known sample geometry from real AFM measurements, with reconstruction quality evaluated as a function of the probe model used. Independent probe characterization (SEM imaging, datasheets, or tip-check samples) provides ground truth for the probe geometry that enters the simulator and against which any probe-estimation method is evaluated.

## 1. Scientific motivation

AFM is widely used because it can resolve nanoscale surface structure without requiring conventional optical imaging. However, the measured topography depends on both the sample and the probe. This creates an interpretation problem:

1. the true sample is generally unknown;
2. the effective probe shape may be imperfect or may change during scanning; and
3. different combinations of sample and probe geometry can produce similar measurements.

The resulting inverse problem is therefore non-unique in general. Some information is never observed because the probe cannot physically enter certain regions. Noise, drift, feedback artifacts, contamination, and uncertainty in the probe shape further reduce recoverability.

The project addresses this challenge in stages. The forward model first makes the image-formation mechanism explicit and testable. It can then generate controlled pairs of true and measured surfaces, which are needed to compare classical reconstruction, blind tip reconstruction, differentiable optimization, and learned restoration methods.

## 2. Scope of the current model

The current software is a geometric, contact-mode height-field simulator. A surface is represented as one height value for each lateral position, and the probe is represented by its lower height envelope. Standard contact-mode imaging is computed as morphological dilation:

`M[i] = max_j { S[i+j] - T[j] }`

where:

- `S` is the true surface;
- `T` is the probe height above its contact reference; and
- `M` is the simulated AFM measurement.

Morphological erosion is also available:

`E[i] = min_j { S[i+j] + T[j] }`

Standard forward imaging uses dilation for every probe geometry, including the concave/inverse-sphere model. Concavity changes the probe kernel and contact reference; it does not automatically change the forward operation to erosion. Erosion is available only when deliberately selected through the inverse controls. It is useful for teaching morphological duality and is the classical starting point for reconstruction when the probe is known. It should not be interpreted as guaranteed recovery of the true surface; it gives the tightest bound consistent with the model, and inaccessible regions remain uncertain.

The simulator is not currently a finite-element model of cantilever dynamics, material deformation, adhesion, friction, or feedback control. It also does not represent overhanging surfaces. These boundaries are deliberate and should be stated clearly whenever results are presented.

## 3. System architecture and workflow

### Browser applications

The project contains two complementary interfaces:

- **2D line-scan simulator:** emphasizes intuition, probe motion, line profiles, local errors, and force-distance visualization.
- **3D raster-scan simulator:** computes an 80 × 80 height map over a 200 × 200 nm field and displays the true and measured meshes using Three.js.

Both interfaces expose probe and surface parameters, noise, transforms, scan controls, error metrics, and export functions. The 3D application additionally supports OBJ import/export and mobile touch interaction.

### Forward-model pipeline

The main computational workflow is:

1. Generate or import a sample surface.
2. Apply deterministic noise and optional geometric transforms.
3. Generate or import the probe kernel.
4. Compute dilation at every scan position, unless erosion is explicitly selected.
5. Calculate the absolute-error map, maximum error, and root-mean-square error.
6. Display the true surface, simulated measurement, and explanatory probe animation.
7. Export data and parameters for independent analysis.

The moving probe is an explanatory visualization layer. The measured height map is computed from the complete discrete grid-based morphological operation, not from the animation. Bilinear surface sampling remains part of the 3D collision-clearance calculation used to place the rendered probe smoothly between grid points; it is not used to calculate the measured height map. Separating these responsibilities prevents visual placement, vertical exaggeration, or animation speed from altering the numerical output.

### Python model

`afm_forward_model.py` mirrors the established 3D browser physics using NumPy arrays and explicit configuration objects. It currently covers the original analytical probe set but has not yet been extended to the newer faceted-pyramid model. It provides a route to:

- automated comparison with browser exports;
- repeatable parameter sweeps;
- regression tests and reference outputs;
- larger synthetic datasets; and
- integration with inverse algorithms and machine-learning frameworks.

### Geometry exchange and validation

The 3D simulator accepts custom probe and surface meshes in Wavefront OBJ format. Imported triangles are rasterized onto the simulator grid:

- the lower envelope becomes the probe kernel; and
- the upper envelope becomes the sample height map.

`validate_tip_geometry.py` independently parses OBJ files, normalizes their coordinate convention, intersects vertical rays with triangles, and compares the resulting lower envelope against the analytical probe equation. It produces machine-readable JSON, sampled CSV data, profile SVGs, and an interactive Three.js comparison.

## 4. Implemented capabilities

The current platform includes the following capabilities.

### Probe models

Across the applications, the analytical library covers conical, spherical, hyperboloidal, flat-punch, concave, asymmetric, double-tip, and combined sphere-cone behaviours. The 3D application additionally includes a faceted-pyramid model with unequal front, back, and side angles and a rounded apex; its default dimensions are informed by a commercial probe specification. Arbitrary height-field-compatible probe meshes can be imported into the 3D application through OBJ.

These models expose characteristic artifacts such as feature broadening, peak rounding, inaccessible valleys, directional displacement, and double-tip ghosting.

### Surface models

The synthetic surface library includes periodic waves, rolling hills, hemispherical particles, pits, trenches, rough surfaces, pyramids, ridge-like structures, chirps, and lattices. Custom OBJ surfaces can also be imported.

The range is intentionally mixed: simple geometries support sanity checks, while multiscale and directional surfaces create more challenging artifact patterns.

### Reproducibility and interoperability

The applications can export:

- true and measured CSV height maps;
- true and measured 16-bit PNG data images;
- JSON parameters and height-decoding metadata;
- current probe and surface meshes as OBJ from the 3D application; and
- grouped exports for complete simulation cases.

The browser interface is responsive, supports light and dark themes, and includes mouse and touch camera controls. It can be deployed as static files and requires no local installation for standard use.

## 5. Key design choices and alternatives considered

### 5.1 Mathematical morphology rather than full contact mechanics

**Choice:** Use height-field dilation for forward imaging, with erosion exposed as an optional inverse operation.

**Reason:** This directly models the dominant geometric tip-convolution effect, follows established AFM image-simulation literature, is interpretable, and is fast enough for interactive use.

**Alternative:** Finite-element or molecular contact simulation.

**Why not at this stage:** Those approaches could represent deformation and force-dependent contact, but they require material properties, boundary conditions, substantially more computation, and a different validation programme. They would obscure the current research question about geometric probe artifacts.

### 5.2 Height fields rather than unrestricted volumetric geometry

**Choice:** Represent the surface and effective probe as single-valued heights over a lateral plane.

**Reason:** AFM topography is naturally stored as a raster height map, and morphological operations are well defined and computationally manageable in this representation.

**Alternative:** Volumetric voxels, signed-distance fields, or general mesh contact.

**Why not at this stage:** General 3D representations can preserve overhangs and sidewalls but greatly increase storage, contact-search complexity, and ambiguity when converting results back to standard AFM images. OBJ support is retained, but meshes are explicitly reduced to their relevant height envelopes.

### 5.3 A browser-first interface rather than a desktop-only tool

**Choice:** Implement the interactive simulator as static HTML, CSS, and JavaScript, with Three.js for 3D rendering.

**Reason:** A zero-install application is easy to demonstrate, distribute, teach with, and deploy. Immediate visual feedback is particularly valuable for explaining a geometric phenomenon.

**Alternative:** A Python desktop interface or specialized scientific software plugin.

**Why not exclusively:** A desktop application would simplify access to scientific libraries but would add environment and packaging requirements. The project instead uses Python as a parallel analysis layer while keeping the interactive experience portable.

### 5.4 Complementary 2D and 3D tools rather than a single interface

**Choice:** Maintain a fast 2D line scan and a spatially complete 3D raster scan.

**Reason:** The 2D view makes the equations and contact geometry easier to understand. The 3D view captures directional, faceted, and non-axisymmetric effects. Each interface serves a different explanatory purpose.

**Alternative:** Use only the 3D simulator.

**Why not:** A 3D view is visually rich but can hide the exact profile relationship and is computationally heavier. The 2D model is also the most defensible first setting for probe profiles extracted from a single side-view SEM image.

### 5.5 An 80 × 80 interactive grid rather than maximum resolution

**Choice:** Use a modest fixed grid for real-time 3D interaction.

**Reason:** The direct morphological computation evaluates a neighbourhood at every pixel. The chosen grid gives responsive parameter changes and animation in a browser while retaining sufficient structure for demonstrations.

**Alternative:** Use 256 × 256 or larger grids throughout.

**Why not in the interactive path:** Runtime increases rapidly for direct neighbourhood searches. Higher-resolution batch simulation belongs in the Python workflow, where vectorization, compiled morphology routines, parallel processing, or GPU acceleration can be introduced and benchmarked.

### 5.6 Analytical probes plus OBJ import rather than either approach alone

**Choice:** Support both parameterized equations and imported meshes.

**Reason:** Analytical probes are reproducible, smooth, and suitable for controlled experiments. OBJ import allows worn, asymmetric, CAD-derived, or SEM-informed approximations to enter the same pipeline.

**Alternative:** Use only ideal equations or only meshes.

**Why not:** Equation-only simulation cannot represent realistic non-ideal probes, while mesh-only simulation makes parameter sweeps and exact reference comparisons more difficult.

### 5.7 A separate Python reference path rather than relying only on browser output

**Choice:** Port the forward physics to Python.

**Reason:** Independent implementations make discrepancies detectable and support automated testing, batch generation, and future inverse modelling.

**Alternative:** Run the JavaScript simulator headlessly for every experiment.

**Why not exclusively:** Headless browser automation would reuse the same implementation and could reproduce the same unnoticed defect. A second implementation provides stronger cross-checking, although both still require external physical validation.

### 5.8 Deterministic synthetic noise rather than uncontrolled randomness

**Choice:** Seed the noise process and record simulation parameters.

**Reason:** Identical inputs can be regenerated for debugging, comparison, and publication.

**Alternative:** Generate fresh random noise on every run.

**Why not:** Unrecorded randomness makes regression failures difficult to diagnose and weakens reproducibility.

### 5.9 Separate numerical measurement from probe placement

**Choice:** Compute the full measured grid independently, then place the rendered probe using collision-clearance sampling.

**Reason:** The animation communicates the process, while the grid operation remains the numerical source of truth. This separation also permits vertical display exaggeration without changing physics.

**Alternative:** Derive output directly from the animated probe position.

**Why not:** Display meshes and animation steps are optimized for perception, not numerical completeness. Coupling them would make output dependent on frame rate, scan animation, and rendering resolution.

### 5.10 Quantitative validation rather than visual agreement alone

**Choice:** Compare mesh envelopes and analytical equations at sampled coordinates, report coverage and error statistics, and use explicit tolerances.

**Reason:** Two 3D shapes can look similar while differing enough to change the forward model. Machine-readable metrics support regression testing and publication evidence.

**Alternative:** Rely on screenshots or manual inspection.

**Why not:** Visual checks are useful for communication but cannot establish numerical consistency.

### 5.11 Self-contained static applications rather than an early framework migration

**Choice:** Keep each simulator self-contained and deployable as a static page during the prototype and demonstration phase.

**Reason:** This preserves zero-install use, minimizes deployment dependencies, and makes a complete demonstration easy to archive or share.

**Alternative:** Split the JavaScript into modules and introduce a package manager, test runner, and build system.

**Why not yet:** A framework migration would not directly answer the current scientific questions. However, the single-file 3D application is now large, so modularization is becoming justified for unit testing, shared browser/Python definitions, and maintainability before dataset-scale or inverse-model development.

## 6. Validation status

Several levels of checking are now available:

- analytical self-checks for probe equations;
- Python syntax checks;
- flat-surface sanity checks, for which a correctly normalized probe should produce zero geometric artifact;
- independent comparison of OBJ lower envelopes with analytical probe profiles; and
- exported CSV, SVG, JSON, and interactive 3D evidence.

Two current mesh-validation examples pass the configured requirements of at least 99% coverage and no sampled absolute error above 0.25 nm:

- The 15° cone achieved 100% sample coverage, an RMSE of approximately 0.117 nm, and a maximum absolute error of approximately 0.224 nm. The residual is consistent with finite polygonal discretization.
- The 15° triangular pyramid achieved 100% sample coverage, an RMSE of approximately `2.36 × 10^-6` nm, and a maximum absolute error below `5 × 10^-6` nm.

These results validate specific exported mesh profiles against their intended equations. They do not yet validate the complete simulator against experimental AFM data.

The agreed next validation steps close the gap between internal consistency and experimental accuracy through two complementary comparisons:

**Probe-geometry validation via independent characterization.** An experienced collaborator will demonstrate the procedure for extracting calibrated raw profile data from SEM images of AFM probes. That work will:

1. establish the SEM image scale and coordinate convention;
2. digitize and preserve raw probe-profile coordinates;
3. document any alignment, baseline, segmentation, or smoothing applied to the data;
4. convert the extracted profile into a simulator-compatible probe representation without losing its physical scale;
5. compare the extracted profile with the simulator's generated or rasterized probe output; and
6. report residual errors and the uncertainty introduced by image resolution and profile extraction.

**Forward-model validation against real AFM measurements.** Real AFM data is available as Bruker `.spm` files, which preserve calibrated height values, scan parameters, and instrument metadata. By scanning a known calibration sample with a characterized probe, the full forward-model prediction loop can be tested: the simulator takes the known sample geometry and independently characterized probe as inputs, computes the predicted measurement, and the result is compared directly against the real `.spm` height map. Discrepancies reveal modelling gaps that internal self-tests cannot detect.

Together, these two comparisons validate the simulator at the probe-geometry level (SEM profile match) and at the image-formation level (predicted vs real measurement). The `.spm` format is preferred over generic `.tiff` export because it retains the instrument's native calibration, scan metadata, and data-processing history, reducing the risk of misaligned scales or undocumented transformations.

After these experimental comparisons, the supporting validation work for the October manuscript is limited to the highest-value checks that fit a 10-hour weekly budget:

1. browser-to-Python golden tests for representative probe/surface combinations;
2. CSV and 16-bit PNG round-trip tests where needed for claimed exports;
3. deliberate negative controls showing that incorrect geometry fails; and
4. targeted resolution-sensitivity checks rather than exhaustive convergence studies.

Broader automated regression suites and full experimental AFM/reference-sample validation remain important, but they are deferred if they threaten the 1 October draft or 31 October submission-ready deadlines.

## 7. Current limitations and research risks

- A height-field model cannot represent overhangs or multiple heights at one lateral position.
- Geometric convolution does not include deformation, adhesion, friction, cantilever response, drift, line flattening, or feedback dynamics.
- The direct grid algorithm trades resolution for browser responsiveness.
- A single SEM side view constrains only a visible 2D profile; it does not uniquely determine a full 3D probe. Additional characterization (tip-check samples, multiple SEM views, or manufacturer specifications) may be needed for full 3D probe modelling.
- An OBJ mesh may contain a visually complete object whose relevant AFM height envelope is sparse or poorly sampled.
- Imported OBJ probes are automatically centered, shifted to a contact reference, and rescaled when their lateral span falls outside expected bounds. This is convenient for demonstration but can alter absolute nanometre calibration unless scale is checked explicitly.
- The interactive 3D model uses a fixed lateral probe-kernel footprint of approximately 28 nm; geometry outside that footprint does not contribute to the calculation.
- Synthetic validation establishes internal consistency, not experimental accuracy.
- The inverse problem is non-identifiable in some regions, so a visually plausible reconstruction may still be unsupported by the measurement.
- The newest probe features are not yet fully synchronized across the 2D interface, 3D interface, Python model, and user guide.
- The batch dataset generator described in the existing user guide is planned but is not present in the current repository.
- Generated validation reports are currently excluded from version control, so publication-grade evidence will require a deliberate archival policy or compact tracked reference set.
- The worn and split-apex demonstration meshes do not have exact analytical reference profiles and therefore require different validation criteria from the current equation-to-mesh checks.
- Classical morphology, synthetic AFM generation, blind tip reconstruction, and ML restoration are already established. Publication claims must therefore focus on a specific evaluated contribution rather than presenting these components as new in isolation.

## 8. Publication and dissemination plan

The publication strategy is staged and evidence-led. The near-term commitment for this semester is a forward-model manuscript with a complete draft circulated by **1 October 2026** and a submission-ready package by **31 October 2026**. Venue selection and authorship details should be confirmed collectively during the first week of drafting.

### Stage A: Software and demonstration output (ongoing support for Stage B)

Prepare a versioned software snapshot with documented coordinate conventions, example geometries, and reproducible validation commands. Archive the SEM-profile comparison evidence, `.spm`-based forward-model comparisons, and mesh-validation reports used in the paper. Posters, demonstrations, and research meetings can continue to use the simulator, but feature expansion is deferred until after the manuscript is submission-ready.

### Stage B: Forward-model paper (complete draft 1 October; submission-ready 31 October)

Freeze the forward-paper scope around the validated forward simulator, experimental validation against real AFM data, probe-geometry comparison, current mesh validation, exports, limitations, and reproducibility. A methods or software paper is credible if the contribution is framed around the integrated and validated platform rather than claiming novelty for morphological dilation itself.

A submission-ready package should include:

- forward-model prediction vs real `.spm` measurement on a known calibration sample, with residual analysis;
- SEM raw-profile extraction, calibration, and residual comparison against the simulator's probe representation;
- current analytical/mesh geometry validation (for example cone and triangular pyramid);
- the highest-value browser-to-Python agreement cases that fit the October budget;
- targeted negative-control or resolution-sensitivity checks where feasible;
- openly documented data, commands, and evaluation protocols.

If work slips, protect the October deadlines by moving additional tip types, broad convergence studies, large dataset generation, major UI refactoring, and software restructuring to post-submission work.

The likely contribution is a transparent bridge between interactive education, arbitrary probe geometry, reproducible forward simulation validated against real measurements, and a path to later dataset generation and inverse benchmarking. Venue selection should follow a brief review of expectations for software, microscopy methods, and engineering-education publications.

### Stage C: Main inverse-problem study (foundation by 27 November; full study next semester)

The longer-term research opportunity is a controlled study of inverse reconstruction grounded in real AFM data. The availability of real `.spm` measurements from known calibration samples, combined with independent probe characterization (SEM, datasheets, or tip-check samples), enables a direct experimental benchmark that most prior work has not attempted. The candidate research question is:

> Can inverse AFM reconstruction methods recover the known geometry of a calibration sample from real AFM measurements, and how does reconstruction quality depend on the probe model used—ideal, parameterized non-ideal, or independently characterized?

This question is grounded in real data rather than a purely synthetic domain-shift study, making the results directly relevant to practitioners. Before committing to specific novelty claims, the focused literature brief should be expanded into a systematic search with database queries and forward/backward citation chaining. This semester ends with a dataset schema, baselines, and a pilot protocol—not a completed inverse paper.

### Stage D: Extended real-data validation (after the semester foundation)

With additional calibration samples and probes, the benchmark can be expanded to evaluate how reconstruction methods behave across different probe conditions (new, worn, contaminated) and sample geometries. Ideally, SEM observations of the probe before and after AFM scanning would provide independent evidence of probe wear and allow pre/post comparisons. This would support a follow-on microscopy application paper and strengthen claims about practical reconstruction value beyond simulation-only evaluation.

Authorship, contribution statements, software licensing, data release, and target venues should be agreed collectively before manuscript preparation.

## 9. Future work: solving the inverse problem

The inverse programme should progress from constrained and interpretable baselines toward more flexible models.

### Phase 1: Establish recoverability and classical baselines

1. Implement known-probe erosion as a reconstruction bound.
2. Add certainty maps to identify locations where reconstruction is supported by a unique-contact condition.
3. Integrate or reproduce a classical blind tip reconstruction baseline.
4. Report surface error and probe error separately; do not reduce evaluation to visual quality.
5. Include naive baselines such as using the measured image unchanged.

This phase establishes what information the forward model permits before introducing machine learning.

### Phase 2: Build a controlled dataset

Generate paired examples containing:

- true surface (known calibration geometry or synthetic ground truth);
- measured AFM surface (simulated via the forward model and, where available, real `.spm` measurements);
- probe kernel (analytical, parameterized, or independently characterized);
- reconstruction or certainty mask;
- noise and artifact settings; and
- complete metadata describing units, grid spacing, transforms, and random seeds.

Probe data should be divided into distinct groups:

- ideal analytical probes;
- parameterized faceted probes;
- synthetic worn or multi-apex probes; and
- independently characterized probes (from SEM profiles, datasheets, or tip-check samples).

Entire probe geometries must be held out during testing. Randomly splitting images made by the same probes would leak probe-specific information and underestimate domain shift. Where real `.spm` data is available for known samples, it should be included as an additional evaluation tier alongside synthetic forward-model pairs.

### Phase 3: Evaluate reconstruction quality across probe models

Evaluate inverse methods by comparing recovered surface geometry against known ground truth. Vary the probe model provided to each method—ideal, parameterized non-ideal, or independently characterized—and measure how reconstruction quality depends on probe-model fidelity. Stratify results by probe condition, surface class, noise, and feature scale.

Where real `.spm` measurements of known calibration samples are available, these provide the strongest evaluation tier: success means recovering the known sample geometry from a real measurement, not just from a synthetic forward-model output.

Recommended metrics include height RMSE and mean absolute error, structural similarity, feature width/height error, spectral error, uncertainty calibration, and probe-kernel error when the probe is estimated.

### Phase 4: Test mitigation methods

Compare, under the same splits:

- classical erosion with a known probe;
- classical blind tip reconstruction;
- differentiable blind tip reconstruction with regularization;
- a supervised restoration baseline;
- ideal-only training;
- ideal plus non-ideal training;
- domain randomization over probe geometry; and
- physics-informed or unrolled models that embed dilation/erosion in the network.

The main ablation should isolate whether using a more accurate probe model improves reconstruction of known sample geometry, not merely whether a larger model improves average error on synthetic data.

### Phase 5: Add uncertainty and joint estimation

Because multiple surface/probe pairs can explain one measurement, the model should eventually estimate uncertainty or a set of feasible solutions rather than only one deterministic image. Promising directions include:

- confidence or certainty maps;
- Bayesian probe and surface estimation;
- ensembles or probabilistic decoders;
- self-consistency through the differentiable forward model; and
- joint estimation from multiple scans, scan directions, or tilted views.

Multi-view or multi-orientation data are especially important for resolving geometry that a single height map or SEM view cannot identify.

### Phase 6: Closed-loop experimental validation

Apply the selected method to real `.spm` measurements and verify results against independently known geometry. The validation loop is:

1. obtain a real AFM measurement of a known calibration structure;
2. characterize the probe independently (SEM, datasheet, or tip-check);
3. reconstruct the surface and, where applicable, the probe from the real measurement;
4. compare the reconstructed surface against the known calibration geometry;
5. forward-simulate the reconstruction to check self-consistency with the measurement; and
6. report discrepancies and uncertainty rather than selecting only visually successful cases.

This loop has been central to the project from the outset: the forward model is validated against real data, and inverse methods are evaluated by their ability to recover known geometry from real data.

## 10. Semester research timeline (8 September – 27 November 2026)

Capacity is approximately **115 hours** over this window at **10 hours per week**: about 8 hours of focused research and 2 hours for meetings, review, documentation, and contingency. Inverse-model training, a large dataset, major UI refactoring, and nonessential simulator features are deferred until after the forward-model paper is submission-ready.

### 10.1 8 September – 1 October: complete manuscript draft

| Window | Hours | Work | Exit criteria |
|--------|------:|------|---------------|
| 8–13 Sep | 10 | Agree on paper contribution and claims; select a likely venue/template; build the manuscript outline and figure/evidence checklist; assign collaborative review responsibilities. | Outline, claim list, and evidence checklist approved. |
| 14–20 Sep | 10 | Learn the SEM raw-data extraction workflow; extract and calibrate one probe profile; compare it with the simulator representation; begin `.spm`-based forward-model comparison on a known calibration sample; preserve raw/processed data and uncertainty notes. | Calibrated SEM profile comparison and initial `.spm` forward-model comparison with residuals and uncertainty notes. |
| 21–27 Sep | 10 | Rerun current cone/pyramid validation; prepare core figures including simulator-vs-real comparison; write Methods, Results, validation, and limitations. | Core figures and Methods/Results/limitations draft complete. |
| 28 Sep – 1 Oct | 5–7 | Assemble Introduction, related work, discussion, conclusion, references, and collaborative contribution wording. | **Complete draft circulated by 1 October.** |

### 10.2 2–31 October: review and submission package

| Window | Hours | Work | Exit criteria |
|--------|------:|------|---------------|
| 2–4 Oct | 3–5 | Collect and triage feedback; lock the revision list and claims. | Revision list frozen. |
| 5–11 Oct | 10 | Address scientific feedback; add the highest-value browser-to-Python golden cases; reconcile code/documentation inconsistencies. | Priority scientific revisions and golden cases complete. |
| 12–18 Oct | 10 | Complete targeted negative-control and resolution-sensitivity checks; finalize SEM-comparison and `.spm` forward-model comparison metrics, uncertainty, and limitations. | **Validation and uncertainty package frozen by 18 October.** |
| 19–25 Oct | 10 | Complete collaborative revision, reference verification, figure polishing, captions, reproducibility statement, and supplementary material. | Near-final manuscript and supplements ready for last review. |
| 26–31 Oct | 8–10 | Final technical and editorial review, venue formatting, author/contribution approval, and archive of the evidence package. | **Manuscript submission-ready by 31 October.** |

**Deadline-protection rule.** If work slips, protect the paper deadline by moving additional tip types, broad convergence studies, and software restructuring to post-submission work.

### 10.3 1–27 November: inverse-problem foundation

| Window | Hours | Work | Exit criteria |
|--------|------:|------|---------------|
| 1–8 Nov | 10 | Extend the literature matrix with inverse methods evaluated on real data; refine the publication question around recovering known sample geometry from real `.spm` measurements as a function of probe-model fidelity. | Updated literature matrix and narrowed candidate claim. |
| 9–15 Nov | 10 | Define a versioned dataset/metadata schema incorporating real `.spm` data and independently characterized probes; synchronize the Python model with priority probe types; produce a small reproducible pilot set. | **Small versioned dataset prototype with real-data tier by 15 November.** |
| 16–22 Nov | 10 | Establish known-probe erosion and one published reconstruction baseline on both synthetic and real data; define surface, probe, feature, and uncertainty metrics. | Executable baseline and metric definitions. |
| 23–27 Nov | 6–8 | Run a pilot reconstruction on real `.spm` data from a known sample if inputs are ready; otherwise finalize the executable benchmark protocol, risk register, data requirements, and next-semester experiment plan. | **Inverse benchmark specification plus baseline/pilot evidence by 27 November.** |

Solving the full inverse problem is explicitly next-semester work.

### 10.4 Semester milestone checklist

- **1 October:** complete forward-model paper draft circulated.
- **18 October:** validation and uncertainty package frozen.
- **31 October:** forward-model manuscript submission-ready.
- **15 November:** small versioned dataset prototype with real-data tier and probe-level splits.
- **27 November:** inverse benchmark specification plus baseline/pilot evidence; full inverse study deferred.

## 11. Criteria for project success

This semester will be successful if the forward-model manuscript is complete by 1 October, submission-ready by 31 October, and the inverse programme ends on 27 November with an executable benchmark foundation rather than an unfinished paper claim.

The forward-model phase will be successful when the numerical outputs are reproducible, independently testable, consistent across browser and Python implementations for the cases claimed in the paper, validated against real AFM measurements from known calibration samples, and bounded by clearly stated assumptions.

The inverse phase will be successful if it demonstrates—positively or negatively—whether inverse methods can recover known calibration-sample geometry from real AFM measurements, how reconstruction quality depends on probe-model fidelity, and whether independently characterized probes improve performance over idealized models. A negative result would still be valuable if the benchmark is rigorous and reveals where idealized probe models are sufficient.

The broader project succeeds by making AFM probe artifacts easier to understand while producing a defensible path from visualization to quantitative reconstruction research.

## Selected references

1. J. S. Villarrubia, “Algorithms for Scanned Probe Microscope Image Simulation, Surface Reconstruction, and Tip Estimation,” *Journal of Research of the National Institute of Standards and Technology*, 102, 425–454, 1997. <https://doi.org/10.6028/jres.102.030>
2. Y. Matsunaga et al., “End-to-end differentiable blind tip reconstruction for noisy atomic force microscopy images,” *Scientific Reports*, 13, 129, 2023. <https://doi.org/10.1038/s41598-022-27057-2>
3. V. Kocur et al., “Correction of AFM data artifacts using a convolutional neural network trained with synthetically generated data,” *Ultramicroscopy*, 246, 113666, 2023. <https://doi.org/10.1016/j.ultramic.2022.113666>
4. L. K. S. Bonagiri et al., “Precise Surface Profiling at the Nanoscale Enabled by Deep Learning,” *Nano Letters*, 2024. <https://doi.org/10.1021/acs.nanolett.3c04712>
5. S. Chen et al., “Multi-view neural 3D reconstruction of micro- and nanostructures with atomic force microscopy,” *Communications Engineering*, 3, 131, 2024. <https://doi.org/10.1038/s44172-024-00270-9>
6. D. Nečas and P. Klapetek, “Gwyddion: an open-source software for SPM data analysis,” *Central European Journal of Physics*, 10, 181–188, 2012.

