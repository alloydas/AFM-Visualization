# LAS Honors Project Proposal — Copy-Paste Draft

**Course:** HONS 4900H  
**Student:** Avish Cooray  
**Project advisor:** Dr. Soumik Sarkar, Mechanical Engineering  
**Intended term:** Spring 2027  

Use this document to fill the LAS Honors Project Proposal Form. Text under each heading is written for a non-specialist faculty committee. Technical terms are defined in ordinary language.

Before pasting:

- Change **"Will this project also contribute to another course or program requirement"** from Yes to **No**, unless this same project also satisfies a second academic requirement (for example a Computer Science or Data Science independent study in addition to honors). HONS 4900H is the honors project itself, not "another" requirement.
- Confirm **3 credits** and the exact course number with Dr. Sarkar and LAS (`lashnrs@iastate.edu`). LAS requires `4900H` for 1 or more credits; the registrar entry may be `HON 4900A`, `COM S 4900H`, or `DS 4900H`.
- Official catalog title for **MATH 2660** is *Elementary Differential Equations*, not Calculus III. Calculus III is **MATH 2650**. List both if you took both.
- Leave question 7 blank. Check **No** for a joint project. Check **"This project does not require outside approval or training."**
- Send Dr. Sarkar the *Guidelines for Project Advisors* with the form. Deadline is **September 21**.

---

## Student and project header

**Name:** Avish Cooray  
**ISU E-Mail:** avish@iastate.edu  
**Phone Number:** (515) 815-0758  
**Major(s):** Data Science and Computer Science  
**Minor(s):** Artificial Intelligence  
**Academic Advisor(s):** Henry Crowley  

**Title of Project:** Recovering Surfaces and Probes from Atomic Force Microscopy Measurements: The Role of Non-Ideal Probe Geometry  

**Type of Project:** Research  

**Project Advisor and Department:** Dr. Soumik Sarkar [Mechanical Engineering]  

**Will this project also contribute to another course or program requirement:** No  

If the Yes box must be kept, use the following explanation and then stop using it if LAS says the correct answer is No:

> This work is being developed in Dr. Sarkar's laboratory and is also intended as a scholarly manuscript. The honors project is not the laboratory's current forward-model paper, which documents a simulator that predicts AFM images from known surfaces and probes. My honors project addresses the inverse problem in both complementary directions: recovering unknown surface geometry when the probe is treated as known, and recovering unknown probe geometry when the true surface is known from a calibration sample. That study requires its own methods, experiments, evaluation, and written contribution.

---

## 1. Research question(s)

Atomic force microscopy (AFM) maps the shape of a surface at the nanoscale by scanning a very small physical probe, often called a tip, across a sample. Because the probe has a finite shape, the recorded image is not a direct copy of the true surface. Features can appear broader or rounder than they are, and narrow gaps may never be reached. The **forward problem** asks what image a known surface and known probe would produce. The **inverse problem** starts from measured data and tries to recover the unknown geometry that produced it.

The inverse problem has two complementary directions. This project addresses both:

1. **Surface reconstruction (known probe, unknown surface).** Given a measured AFM image and a probe model, can reconstruction methods recover the known geometry of a calibration sample? How does reconstruction quality depend on the probe model used: an idealized shape, a parameterized non-ideal shape, or a probe characterized independently, for example from electron-microscope images or manufacturer specifications?
2. **Probe reconstruction (known surface, unknown probe).** Given a measured AFM image and a known true surface, can the same inverse framework recover the probe geometry? How close is that estimate to an independent description of the same probe, such as a manufacturer specification or an electron-microscope profile?

A central question across both directions is whether using realistic, non-ideal probe geometries during reconstruction or training improves results compared to relying on idealized probe shapes. I hypothesize that each recovery is more accurate when the probe model is closer to the physical probe, and that some regions will remain unrecoverable because the probe and sample never made a unique contact. A negative result would still be useful if it shows where idealized probe models are sufficient and where they fail.

---

## 2. Short description of the project

*(Shorter version — 100 words)*

Atomic force microscopy (AFM) maps nanoscale surface shape by scanning a physical probe across a sample. Because the probe has a finite shape, the recorded image mixes the true surface with the probe. This project studies two inverse problems: recovering the true surface when the probe is treated as known, and recovering the probe when the true surface is known. Existing reconstruction methods typically assume idealized probe shapes. Using a previously developed forward simulator and measurements of known calibration samples, I will test whether training or reconstructing with realistic, non-ideal probe geometries improves recovery compared to idealized models.

---

## 3. Difference from previous work / unique contribution

The geometric relationship between AFM images, surfaces, and probes is well established. Villarrubia (1997) showed that contact-mode imaging can be modeled as a morphological operation, and that the inverse can be run in either direction: a known probe can be used to bound the surface, and a known calibration surface can be used to bound the probe. Software tools such as Gwyddion already implement these classical methods.

More recently, machine-learning approaches have been applied to AFM reconstruction. Kocur et al. (2023) trained a neural network on synthetic AFM data to correct imaging artifacts, and Bonagiri et al. (2024) trained an encoder-decoder network to recover nanoscale surface profiles. Matsunaga et al. (2023) introduced a differentiable version of the classical probe-estimation method that is more robust to noise. However, these learning-based methods generate their training data using idealized probe shapes: perfect cones, spheres, or simple parameterized models. Real probes are never perfectly ideal — they are asymmetric, faceted, worn, or contaminated, and the imaging artifacts they produce differ qualitatively from those of ideal probes.

The specific question this project tests is whether reconstruction methods trained or configured with non-ideal, independently characterized probe geometries perform better than the same methods trained with idealized probes alone. My laboratory has already developed a forward simulator capable of generating synthetic AFM data with a library of non-ideal probe models — including asymmetric, faceted, and multi-apex geometries — as well as real AFM measurements of known calibration samples. That infrastructure makes a controlled comparison possible.

This work is part of a larger laboratory program on AFM image formation and reconstruction in the Translational AI Center, advised by Dr. Soumik Sarkar, with technical guidance from Dr. Adarsh Krishnamurthy and day-to-day mentorship from graduate researcher Alloy Das.

**We** (the laboratory) are developing a validated forward model: software that predicts the AFM image that should be produced by a known surface and probe, including interactive simulators, a Python implementation, and comparison with real AFM data. That forward-model study is a separate paper and is not my honors project. **We** also provide experimental AFM measurements, calibration samples, and, where available, independent probe characterization.

**I** will carry out the inverse study. Specifically, I will:

- establish classical baselines for both inverse directions: surface reconstruction with a known probe, and probe reconstruction with a known calibration surface;
- generate synthetic training data using the forward model with both idealized and non-ideal probe geometries, holding out entire probe shapes during testing;
- train a learning-based reconstruction model under Dr. Sarkar's supervision on idealized probes only, and separately on a mixture that includes non-ideal probes;
- compare the two training regimes on held-out non-ideal probes, measuring whether non-ideal training improves surface and probe recovery;
- apply the methods to real measurements of known calibration samples where the laboratory's data allow;
- report where each recovery succeeds, where information is fundamentally lost, and whether the non-ideal training regime helps;
- write the inverse-study results as my honors contribution.

The forward simulator is infrastructure for this project; it is not itself the honors contribution. The novel experiment is the controlled comparison of reconstruction quality as a function of probe-model realism during training.

---

## 4. Bibliography

Villarrubia, J. S. (1997). Algorithms for scanned probe microscope image simulation, surface reconstruction, and tip estimation. *Journal of Research of the National Institute of Standards and Technology, 102*(4), 425–454. https://doi.org/10.6028/jres.102.030

Matsunaga, Y., Hashimoto, S., Ngo, K. X., Ueno, A., Kodera, N., & Ando, T. (2023). End-to-end differentiable blind tip reconstruction for noisy atomic force microscopy images. *Scientific Reports, 13*, 129. https://doi.org/10.1038/s41598-022-27057-2

Kocur, V., Hegrová, V., Piastek, J., Fiala, P., Neuman, J., & Bábek, O. (2023). Correction of AFM data artifacts using a convolutional neural network trained with synthetically generated data. *Ultramicroscopy, 246*, 113666. https://doi.org/10.1016/j.ultramic.2022.113666

Bonagiri, L. K. S., Wang, Z., Zhou, S., & Zhang, Y. (2024). Precise surface profiling at the nanoscale enabled by deep learning. *Nano Letters*. https://doi.org/10.1021/acs.nanolett.3c04712

Chen, S., et al. (2024). Multi-view neural 3D reconstruction of micro- and nanostructures with atomic force microscopy. *Communications Engineering, 3*, 131. https://doi.org/10.1038/s44172-024-00270-9

Nečas, D., & Klapetek, P. (2012). Gwyddion: An open-source software for SPM data analysis. *Central European Journal of Physics, 10*(1), 181–188.

The form asks for 4–6 sources. All six above are relevant. If a shorter list is required, keep Villarrubia (1997), Matsunaga et al. (2023), Kocur et al. (2023), Bonagiri et al. (2024), and Nečas and Klapetek (2012).

---

## 5. Coursework and experiences

**Coursework completed**

- DS 2010: Introduction to Data Science
- Mathematics through the calculus sequence and MATH 2660: Elementary Differential Equations

**Coursework to be completed before Spring 2027**

- COMS 2280: Introduction to Data Structures
- COMS 2300: Discrete Computational Structures
- STAT 2010: Introduction to Statistical Concepts and Methods
- STAT 3030: Probability and Statistics for Computer Science
- MATH 2070: Matrices and Linear Algebra

Together these courses provide the programming, discrete mathematics, probability, and linear algebra needed to implement reconstruction methods, train models, and evaluate error.

**Research experience**

- *First-Year Honors Mentor Project (Spring 2026, NexGen PPT Lab).* I built a machine-learning pipeline that predicts smoke particle profiles from combustion experiments. The project involved designing a six-module Python data pipeline, engineering features from multi-instrument experimental data, training a Gaussian process surrogate model that improved prediction accuracy from R² 0.026 to 0.759, and deploying the full pipeline to Iowa State's Nova HPC cluster. This experience gave me hands-on practice with the full ML workflow — data cleaning, feature engineering, model training, evaluation metrics, and computational deployment — which directly transfers to training reconstruction models on synthetic AFM data in the proposed project.
- *Summer 2026 undergraduate research (Translational AI Center).* I built a browser-based AFM forward simulator, implemented the corresponding Python model, and presented a research poster on visualizing probe-induced imaging artifacts. That work produced the forward simulator and probe library that this project will use to generate training data for the inverse study; it is background infrastructure, not the honors project itself.

---

## 6. Procedures and methods

All technical terms are defined here in ordinary language.

**Data.** I will use two kinds of data. First, synthetic triples created by the laboratory's forward simulator: a known true surface, a known probe, and the AFM image that probe would record. These triples will be generated with both idealized probes (such as perfect cones or spheres) and non-ideal probes (asymmetric, faceted, or worn shapes drawn from the simulator's probe library). Second, real AFM measurements of known calibration samples (fabricated surfaces whose dimensions are independently specified), together with an independently characterized probe when available. Real files will be used in the instrument's native format so that height calibration and scan settings are preserved.

**Step 1 — Classical baselines for both inverse directions.** Surface reconstruction and probe reconstruction are mathematically related operations that share nearly all the same code. For surface reconstruction, the classical method computes the tallest surface consistent with the measurement and the known probe. For probe reconstruction, the same framework computes the bluntest probe consistent with the measurement and the known calibration surface. I will implement both, verify them on synthetic data, and use them as the baselines that any more flexible method must beat.

**Step 2 — Probe-model comparison with classical methods.** I will run surface reconstruction using different probe models for the same measurement: an idealized shape, a more realistic parameterized shape, and an independently characterized probe. This tests how much the classical result improves when the probe description becomes more accurate, and establishes whether the probe model matters before introducing machine learning.

**Step 3 — Learning-based reconstruction with non-ideal training data.** Under Dr. Sarkar's supervision, I will train a learning-based reconstruction model on synthetic data. The key experiment is a controlled comparison: one model trained on data generated with idealized probes only, and one trained on data that includes non-ideal probe geometries. Both models are evaluated on the same held-out set of non-ideal probes. This tests whether exposure to realistic probe shapes during training improves reconstruction quality.

**Evaluation.** For surface reconstruction, I will measure height error, feature-size error, and maps of which pixels are actually supported by the data. For probe reconstruction, I will measure error in the recovered probe kernel and, where possible, agreement with independent probe characterization. Entire probe shapes will be held out during testing so that a method is not evaluated on probes it was trained or developed on. I will include negative controls, such as reconstructing the surface with the wrong probe or reconstructing the probe from a surface that is too smooth to constrain the tip.

This project does not involve human participants.

---

## 7. Creative/innovative projects only

Not applicable.

---

## 8. Proposed project schedule

I will complete the honors project in Spring 2027, the term of the Honors Poster Presentation, under Dr. Sarkar's supervision. Weekly meetings with my project advisor and laboratory mentors will continue throughout.

| When | What I will complete |
|---|---|
| October–December 2026 (after proposal approval) | Literature review of inverse AFM methods and non-ideal probe studies. Freeze research questions, dataset plan, and evaluation metrics. Confirm which calibration measurements and probe characterizations will be available. Begin implementing the classical baselines. |
| January 2027 | Complete classical baselines for both directions on synthetic data. Run the probe-model comparison (Step 2): test surface reconstruction with idealized, parameterized, and independently characterized probes to quantify how much the probe model matters. |
| February 2027 | Generate training datasets with idealized and non-ideal probe geometries. Train the learning-based reconstruction model under Dr. Sarkar's supervision. Run the controlled comparison: ideal-only training versus non-ideal training, evaluated on held-out probes. |
| March 2027 | Apply methods to real calibration measurements where the laboratory's data allow. Complete residual and uncertainty analysis. Draft results, limitations, and figures. |
| Early April 2027 | Finish the written honors report. Prepare the Honors poster. Circulate a draft to Dr. Sarkar for revision. |
| Late April / early May 2027 | Present at the University Honors Poster Presentation. Archive code, data descriptors, and the written report. |

If the learning-based comparison proves infeasible within the semester, the classical baselines and the probe-model comparison (Steps 1–2) are already a complete scholarly project with clear results. The machine-learning experiment would then be reported as planned future work.

---

## 9. Additional sharing of results

In addition to the Honors Poster Presentation, I will:

- submit a written research report to my project advisor;
- present results to the laboratory research group;
- prepare the inverse-study findings for a scientific manuscript, distinct from the laboratory's forward-model paper;
- deposit a project abstract in the Honors digital repository as required.

---

## 10. Credits and course number

**Number of Credits:** 3  

**Course number:** HONS 4900H  

**Term of Honors poster presentation:** Spring 2027  

---

## 11. Completing the project with another student?

**No**

---

## Responsible research / IRB

**This project does not require outside approval or training.**

The project uses AFM measurements of fabricated calibration samples and computational analysis of those measurements. It does not involve human participants.

---

## Signature checklist

- [ ] I have read the Project Proposal Instructions and agree to the requirements for the completion of this project. *(Student)*
- [ ] I have read the student's proposal and the Guidelines for Project Advisors and am willing to mentor the student in the completion of this project. *(Project advisor: Dr. Soumik Sarkar)*
- [ ] I am aware of my advisee's intention to complete an Honor Project as described in the proposal. *(Academic advisor: Henry Crowley)*

Email the completed form to `lashnrs@iastate.edu` by September 21.
