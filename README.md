# Cortical timescale hierarchy model

Code and derived data supporting the manuscript on predicting intact-versus-sentence-scrambled cortical response variance from a scale-covariant / scale-free filter-bank model.

## What this repository reproduces

The repository reproduces the reported analysis from the **derived participant-level PC1 variances**:

- intact and sentence-scrambled PC1 variance for 24 participants;
- intact-to-scrambled variance ratios;
- geometric-mean variance ratio and bootstrap 95% confidence interval;
- permutation, binomial sign, and paired t-tests;
- the model prediction at the literature-derived sentence-level timescale;
- sensitivity to the input spectral exponent and scrambled-condition maximum integration timescale;
- the model parameter-space figure.

The raw fMRI dataset is not redistributed here.

## Repository structure

```text
.
├── README.md
├── code/
│   ├── run_results.m
│   └── make_parameter_space_figure.m
├── data/
│   └── pc1_variances.csv
└── figures/
```

`run_results.m` reproduces the numerical results and printed tables.

`make_parameter_space_figure.m` computes the model prediction across parameter space and saves PNG, PDF, and SVG versions of the figure to `figures/`.

`pc1_variances.csv` contains the derived participant-level PC1 variances used in the statistical comparison.

## Requirements

- MATLAB R2019b or later
- Statistics and Machine Learning Toolbox

## Running the analysis

Clone or download the repository, open MATLAB, and run:

```matlab
run('code/run_results.m')
run('code/make_parameter_space_figure.m')
```

The scripts use repository-relative paths; no Desktop-specific folder structure is required.

## Public fMRI dataset

The empirical data originate from the Narratives collection, specifically the *It's Not the Fall That Gets You* dataset:

- Nastase et al. (2021), *The Narratives fMRI dataset*, Scientific Data 8, 250.
- OpenNeuro dataset: ds002345.

## Important reproducibility scope

The included CSV contains the **derived PC1 variance values** used in the manuscript. The supplied code therefore reproduces the statistical tests, model predictions, sensitivity analyses, and parameter-space figure from those values.

It does **not yet reproduce the complete path from raw fMRI volumes to the PC1 variance values**. To support a claim that all results can be reproduced directly from the public raw dataset and accompanying code, the repository should additionally contain the preprocessing / ISC-mask / PCA extraction code that generated `data/pc1_variances.csv`.

## Core model parameters

The scripts use:

- minimum integration timescale: 0.05 s;
- maximum intact integration timescale: 45 s;
- reference maximum sentence-scrambled integration timescale: 7.5 s;
- reference input spectral exponent: alpha = 0;
- frequency band: 0.0078–0.1 Hz;
- bootstrap resamples: 10,000;
- permutations: 100,000;
- random seed: 42.

The joint parameter-space analysis examines alpha from 0 to 2 and the scrambled-condition maximum integration timescale from 3 to 12 s.

## References used for model parameters / context

- Honey et al. (2012). *Slow cortical dynamics and the accumulation of information over long timescales*. Neuron 76, 423–434.
- Lerner et al. (2011). *Topographic mapping of a hierarchy of temporal receptive windows using a narrated story*. Journal of Neuroscience 31, 2906–2915.
- Nastase et al. (2021). *The Narratives fMRI dataset*. Scientific Data 8, 250.
- Chang et al. (2022). *Information flow across the cortical timescale hierarchy during narrative construction*. PNAS 119, e2209307119.

## License

No license has been assigned in this package. Add the license you want before publishing the repository publicly.
