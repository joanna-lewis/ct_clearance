# Duration of untreated _C. trachomatis_ infection in men

This repository contains code and other supplementary information for evidence synthesis of the duration of untreated genital _C. trachomatis_ infection in men.

* The directory named `search_strategy` contains details of the search strategy used for updating reviews of the evidence, and a supplementary table (`data.pdf`) summarising all the data used, from men and women.

* The directory named `analysis` contains:
	- Various STAN model specification files;
	- Various datasets, saved as RData objects;
	- Jupyter notebooks for running the STAN models and interpreting output, and
	- A directory named `results` containing saved output from the Jupyter notebooks.
	
Three notebooks are provided:

* `results_2comp.ipynb` runs a model in which durations of chlamydia infection have a two-component mixture-of-exponentials distribution.
* `model_comparison.ipynb` compares the fit of the two-component mixture model with the fits of various alternative models.
* `sensitivity_analysis.ipynb` investigates the sensitivity of the model to the chlamydia diagnosis method used, and to the exclusion of each of the studies in turn.