This repo is included as a submodule in the following two repos:
1. [e-lipman/XXX](<https://github.com/e-lipman/XXX>), which reproduces the results in
[Explaining Differences in Voting Patterns Across Voting Domains Using Hierarchical Bayesian Models (Lipman, Moser, & Rodriguez, 2023)](<https://arxiv.org/abs/2312.15049>)
3. [e-lipman/ModularBayesUSHouse](<https://github.com/e-lipman/ModularBayesUSHouse>), which reproduces the results in
[On Data Analysis Pipelines and Modular Bayesian Modeling (Lipman & Rodriguez, 2024)](https://arxiv.org/abs/2402.04461)

## Data:

## R scripts:
*`1_run_models`: Runs a single chain of a given model for a single House.*

- Inputs (Passed by order from command line)

    -`cong`: A congress number in 93,...,113

    -`chain`: A number to identify the chain

    -`stage`: 0=joint model, 1=stage 1 (working posterior), 2=stage 2 (two-step or cut model)

    -`burn`: Number of burn-in iterations

    -`iter`: Number of post burn-in iterations to save

    -`thin`: Factor by which iterations are thinned

    -`folder`: Name for top level output folder (e.x. "stage1")

-Additional imputs for stage 2 of two-step model and cut models

    -`cut`: Indicator for the cut model (cut=1) versus two-step model (cut=0)

`2_postprocess.R `: Postprecesses results for a given model for a single House. Only needed for stage=0 and stage=1

## C scripts:
