This repo is included as a submodule in the following two repos:
1. [e-lipman/XXX](<https://github.com/e-lipman/XXX>), which reproduces the results in
[Explaining Differences in Voting Patterns Across Voting Domains Using Hierarchical Bayesian Models (Lipman, Moser, & Rodriguez, 2023)](<https://arxiv.org/abs/2312.15049>) (obtains and interprets substantive results using the full model)
3. [e-lipman/ModularBayesUSHouse](<https://github.com/e-lipman/ModularBayesUSHouse>), which reproduces the results in
[On Data Analysis Pipelines and Modular Bayesian Modeling (Lipman & Rodriguez, 2024)](https://arxiv.org/abs/2402.04461) (compares the full, cut, and two-step models)

## Data:

The input data for the XXX-th house is contained in the file `Data/inputs_XXX.RDS`. 

This object contains the following elements used in the C++ modeling code:
    
- `y`: Outcome variable. Matrix of votes (1=yay, 0=nay, NA=absent/abstained). `y` has one row per legislator and one column per vote. Columns are sorted so that the procedural votes are grouped tpgether before the final passage votes

- `missing`: Boolean matrix of same dimension as `y` indicating which elements of `y` are `NA`

- `covariates`: A data frame with one row per house member (same order as rows of `y`) containing covariates for each member (used to create `X` matrix of covariates)

- `gam0_maxC`: Index of last procedural vote using a base-0 index

- `demLeaderC'/`repLeaderC`: Indices of Democratic and Republican party leaders, using a base-0 index


## R scripts:
*`1_run_models`: Runs a single chain of a given model for a single House.*

- Arguments (passed by position as command line arguments)
    
    - For full model (stage=0) and stage 1 (stage=1): "cong", "chain", "stage", "burn", "iter", "thin", "folder"
    
    - For stage 2 of two-step model (stage=2, cut=0): "cong", "chain", "stage", "cut", "burn", "iter", "thin", "folder", "infolder", "insuffix","thresh"
    
    - For stage 2 of cut model (stage=2, cut=1): "cong", "chain", "stage", "cut", "steps", "folder", "infolder", "insuffix"

- Description of arguments:

    -`cong`: A congress number in 93,...,113

    -`chain`: A number to identify the chain

    -`stage`: 0=joint model, 1=stage 1 (working posterior), 2=stage 2 (two-step or cut model)

    -`folder`: Name for top level output folder (e.x. "stage1")

    -`burn`: (not used for stage 2 cut) Number of burn-in iterations

    -`iter`: (not used for stage 2 cut) Number of post burn-in iterations to save

    -`thin`: (not used for stage 2 cut) Factor by which iterations are thinned

    -`cut`: (stage 2 only) Indicator for the cut model (cut=1) versus two-step model (cut=0)
  
    -`infolder`: (stage 2 only) Top=level folder from which to lead stage 1 results (e.x. "stage1")

    -`insuffix`: (stage 2 only) Suffix ("_burn_iter_thin") for stage 1 results

    -`thresh`: (stage 2 of two-step only) Posterior probability threshold for determining bridges from stage 1 (e.x. 0.5)

    -`steps`: (stage 2 of cut only) Number of steps in MCMC chain for each iteration of stage 1

`2_postprocess.R `: Postprecesses results for a given model for a single House. Only needed for stage=0 and stage=1

- Arguments (passed by position as command line arguments): "cong", "stage", "folder", "suffix"

## C++ scripts (in `src`):
  
-`run_sampler.cpp`: Top-level C script to run sampler

-`gibbs_updates.cpp`: Gibbs updates for module 1 

-`covaraite_regression.cpp`: Gibbs/MH upodates for module 2

-`helper_funstions.cpp`: Helper functions

-`parameters.hpp`: Defines structure to hold parameters

-`configs.hpp`: Static hyperparameters for the model, especially parameters and settings for hyperpriors

-`helloPG`: This folder contains files for generating Polya-Gamma random variables. These files are taken from 
[jgscott/helloPG](<https://github.com/jgscott/helloPG>) (current as of 2/2024)
