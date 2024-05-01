# Simulation Study of Generalizability Estimators
__BST 258 Final Project (Spring 2024)__
* Sarika Aggarwal
* Keith Barnatchez
* Luke Benz
* Zhu (Lucy) Shen

---
This repository contains code for a series of simulation studies to understand the properties of various generalizability estimators. It is based on the work of [(Colnet et al., 2024)](https://arxiv.org/abs/2011.08047). Descriptions of relevant files are as follows


* __generalizability_sims.R__: Main simulation script
* __estimators.R__: Implementation of several generalizability estimators using both parametric (LM/GLM) and non-parametric (random forest) models.
* __helpers.R__:  Script defining useful helper functions
* __dgp_figures.R__: Script to produce some figures on both the data generating process and the distribution of outcome predictions for regression based models.
* __shell/__: Folder of shell scripts to transfer files to/from the Harvard computing cluster.
* __sim_results/__: Folder of results. File __generalizability.csv__ contains the results of our primary simulation study.
* __figures/__: Folder of figures.

