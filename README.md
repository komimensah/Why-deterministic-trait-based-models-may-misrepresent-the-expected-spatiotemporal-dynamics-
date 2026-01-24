This repository contains the code and figures supporting a trait-based mechanistic framework for modelling tick population growth potential under climate forcing, with explicit propagation of biological uncertainty.

The model extends a deterministic reproduction number (R₀) formulation by embedding uncertainty in key life-history traits—such as fecundity and host-finding time—directly within the nonlinear population growth structure using Monte Carlo sampling. Deterministic and stochastic formulations are evaluated consistently across space and time using identical climate inputs, allowing direct assessment of how uncertainty alters expected R₀ magnitude, seasonal dynamics, and spatial risk patterns.

Repository contents
	•	Code 1: Trait-based R₀ model with stochastic uncertainty propagation
	•	Code 2: Spatial implementation and quarterly aggregation of deterministic vs stochastic R₀
  • Code 3: Main model to generate the raster outputs
  •	Figure 1: Workflow
	•	Figure 2: Temperature–R₀ relationships illustrating deterministic and uncertainty-propagated expectations
	•	Figure 3: Quarterly spatial comparison of deterministic and stochastic R₀ across North Africa

Key features
	•	Explicit treatment of parametric biological uncertainty within a mechanistic framework
	•	Preservation of biologically plausible spatial patterns while refining expected population growth
	•	Identification of critical seasonal windows relevant for surveillance and intervention planning
	•	Designed for decision-support in data-limited settings

The code is fully reproducible and intended to support transparency, reuse, and further methodological development in uncertainty-aware ecological and epidemiological modelling.
