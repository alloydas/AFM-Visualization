# AFM Forward Simulation, Validation, and Inverse Reconstruction Roadmap

**Project status:** Active research and software development  
**Current emphasis:** Validated forward modelling and preparation for inverse reconstruction  


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

The immediate research objective is to establish a reliable and reproducible forward model. The longer-term objective is to study the inverse problem: recovering the surface, the probe, or both from a distorted AFM measurement. A promising publication direction is not simply “using machine learning for AFM deconvolution,” which is already established, but testing whether inverse methods trained on ideal probe shapes remain reliable for held-out, SEM-informed non-ideal probes.

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

The agreed next validation step is to extract calibrated raw profile data from an SEM image, with the extraction procedure demonstrated by an experienced collaborator. That work will:

1. establish the SEM image scale and coordinate convention;
2. digitize and preserve raw probe-profile coordinates;
3. document any alignment, baseline, segmentation, or smoothing applied to the data;
4. convert the extracted profile into a simulator-compatible probe representation without losing its physical scale;
5. compare the extracted SEM profile with the simulator's generated or rasterized probe output; and
6. report residual errors and the uncertainty introduced by image resolution and profile extraction.

This comparison provides direct validation of the simulator's probe-geometry representation. Validating the complete forward image-formation model will additionally require a known or independently characterized sample, an AFM measurement made with the corresponding probe, and comparison with the simulator's predicted measured surface.

After the SEM profile comparison, the supporting validation work is:

1. browser-to-Python golden tests for representative probe/surface combinations;
2. CSV and 16-bit PNG round-trip tests;
3. deliberate negative controls showing that incorrect geometry fails;
4. convergence studies over grid spacing and probe-mesh resolution; and
5. automated regression tests for transforms, boundaries, and imported geometry.

## 7. Current limitations and research risks

- A height-field model cannot represent overhangs or multiple heights at one lateral position.
- Geometric convolution does not include deformation, adhesion, friction, cantilever response, drift, line flattening, or feedback dynamics.
- The direct grid algorithm trades resolution for browser responsiveness.
- A single SEM side view constrains only a visible 2D profile; it does not uniquely determine a full 3D probe.
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

The publication strategy should be staged and evidence-led. No target venue or manuscript timeline is currently fixed in the project files, so the stages below are a proposed plan for collective review rather than an existing commitment.

### Stage A: Software and demonstration output

Prepare a versioned software release with a stable web deployment, installation-free demonstration, user guide, example geometries, documented coordinate conventions, and reproducible validation commands. Archive a citable release if project policy permits. Posters, demonstrations, and research meetings can use the simulator to communicate the problem and gather feedback.

### Stage B: Forward-model and benchmark paper

A methods or software paper is possible if the contribution is framed around the integrated and validated platform rather than claiming novelty for morphological dilation itself. A credible submission would include:

- formal browser-to-Python agreement tests;
- sensitivity and convergence studies;
- validation of analytical and imported probe geometries;
- a reproducible artifact benchmark spanning ideal and non-ideal probes;
- comparison with an established package such as Gwyddion where applicable;
- at least one physical reference-sample case if data are available; and
- openly documented data and evaluation protocols.

The likely contribution would be a transparent bridge between interactive education, arbitrary probe geometry, reproducible forward simulation, and dataset generation. Venue selection should follow a review of expectations for software, microscopy methods, and engineering-education publications.

### Stage C: Main inverse-problem study

The stronger research publication opportunity is a controlled study of probe-geometry domain shift:

> Do inverse AFM methods trained or tuned with idealized probes retain accuracy on measurements formed by held-out SEM-informed non-ideal probes, and does training with a realistic probe distribution improve robustness?

This wording avoids unsupported claims that SEM tip characterization or ML deconvolution is itself new. Before using a novelty claim, the focused literature brief should be expanded into a systematic search with database queries and forward/backward citation chaining.

### Stage D: Experimental extension

If suitable data can be collected, evaluate the best methods using a known reference sample and an independently imaged probe, ideally with SEM observations before and after AFM scanning. This would distinguish simulation consistency from practical reconstruction value and could support a follow-on microscopy application paper.

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

- true surface;
- measured AFM surface;
- probe kernel;
- reconstruction or certainty mask;
- noise and artifact settings; and
- complete metadata describing units, grid spacing, transforms, and random seeds.

Probe data should be divided into distinct groups:

- ideal analytical probes;
- parameterized faceted probes;
- synthetic worn or multi-apex probes; and
- SEM-informed non-ideal profiles or explicitly labelled 3D approximations.

Entire probe geometries must be held out during testing. Randomly splitting images made by the same probes would leak probe-specific information and underestimate domain shift.

### Phase 3: Evaluate domain shift

Train or configure methods using ideal probes, then evaluate without adjustment on held-out non-ideal probes. Stratify results by probe wear, asymmetry, apex multiplicity, noise, surface class, and feature scale.

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

The main ablation should isolate whether realistic probe distributions improve generalization, not merely whether a larger model improves average error.

### Phase 5: Add uncertainty and joint estimation

Because multiple surface/probe pairs can explain one measurement, the model should eventually estimate uncertainty or a set of feasible solutions rather than only one deterministic image. Promising directions include:

- confidence or certainty maps;
- Bayesian probe and surface estimation;
- ensembles or probabilistic decoders;
- self-consistency through the differentiable forward model; and
- joint estimation from multiple scans, scan directions, or tilted views.

Multi-view or multi-orientation data are especially important for resolving geometry that a single height map or SEM view cannot identify.

### Phase 6: Sim-to-real validation

Apply the selected method to experimental data only after synthetic performance and failure modes are understood. Where possible:

1. scan a known calibration structure;
2. characterize the probe independently;
3. reconstruct the surface and, where applicable, the probe;
4. forward-simulate the reconstruction;
5. compare the re-simulated image with the measurement; and
6. report discrepancies and uncertainty rather than selecting only visually successful cases.

## 10. Proposed near-term milestones

### Milestone 1: Forward-model verification

- Learn and document the raw SEM profile-extraction procedure.
- Extract calibrated probe coordinates and retain both raw and processed data.
- Compare the SEM-derived profile against the simulator representation, including residuals and uncertainty.
- Add browser/Python golden cases.
- Add negative controls and resolution-convergence tests.
- Validate exported faceted-pyramid geometry and define mesh-only checks for worn probes.
- Reconcile documentation with the currently implemented probe library.
- Record all validation commands and software versions.

### Milestone 2: Dataset prototype

- Build a batch generator around the Python model.
- Define a versioned metadata schema.
- Produce a small balanced dataset with probe-level train/validation/test splits.
- Verify every sample by re-running the saved parameters.

### Milestone 3: Inverse baseline benchmark

- Implement known-probe erosion and certainty maps.
- Select classical and differentiable blind-tip baselines.
- Establish quantitative surface and probe metrics.
- Publish a benchmark protocol before tuning a new model.

### Milestone 4: SEM-informed pilot

- Extract defensible 2D profiles from several SEM images.
- Record calibration, segmentation, fitting uncertainty, and hidden-face assumptions.
- Test ideal-to-non-ideal domain shift first in 2D.
- Proceed to faceted 3D approximations only if the pilot demonstrates a measurable and scientifically useful gap.

### Milestone 5: Experimental and publication package

- Add import support for real SPM height maps with documented calibration and preprocessing.
- Collect or obtain reference-sample AFM data.
- Complete ablations and uncertainty analysis.
- Release code, dataset documentation, and reproducible evaluation scripts as permitted.
- Prepare the forward-platform and inverse-domain-shift outputs for the most appropriate venues.

## 11. Criteria for project success

The forward-model phase will be successful when the numerical outputs are reproducible, independently testable, consistent across browser and Python implementations, and bounded by clearly stated assumptions.

The inverse phase will be successful if it determines—positively or negatively—whether realistic probe geometry causes a meaningful generalization failure, and whether a controlled mitigation improves performance on genuinely held-out probes. A negative result would still be valuable if the benchmark is rigorous and reveals where idealized models are sufficient.

The broader project succeeds by making AFM probe artifacts easier to understand while producing a defensible path from visualization to quantitative reconstruction research.

## Selected references

1. J. S. Villarrubia, “Algorithms for Scanned Probe Microscope Image Simulation, Surface Reconstruction, and Tip Estimation,” *Journal of Research of the National Institute of Standards and Technology*, 102, 425–454, 1997. <https://doi.org/10.6028/jres.102.030>
2. Y. Matsunaga et al., “End-to-end differentiable blind tip reconstruction for noisy atomic force microscopy images,” *Scientific Reports*, 13, 129, 2023. <https://doi.org/10.1038/s41598-022-27057-2>
3. V. Kocur et al., “Correction of AFM data artifacts using a convolutional neural network trained with synthetically generated data,” *Ultramicroscopy*, 246, 113666, 2023. <https://doi.org/10.1016/j.ultramic.2022.113666>
4. L. K. S. Bonagiri et al., “Precise Surface Profiling at the Nanoscale Enabled by Deep Learning,” *Nano Letters*, 2024. <https://doi.org/10.1021/acs.nanolett.3c04712>
5. S. Chen et al., “Multi-view neural 3D reconstruction of micro- and nanostructures with atomic force microscopy,” *Communications Engineering*, 3, 131, 2024. <https://doi.org/10.1038/s44172-024-00270-9>
6. D. Nečas and P. Klapetek, “Gwyddion: an open-source software for SPM data analysis,” *Central European Journal of Physics*, 10, 181–188, 2012.

