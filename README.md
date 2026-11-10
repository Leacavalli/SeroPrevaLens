<p align="center">
  <img src="assets/logo.svg" alt="SeroPrevaLens logo" width="120">
</p>

<h1 align="center">SeroPrevaLens</h1>
<p align="center"><b>Estimate SeroPrevalence using Bayesian Mixture Models</b></p>

<p align="center">
  <a href="https://leacavalli.github.io/SeroPrevaLens/"><b>Open the tool</b></a>
</p>

SeroPrevaLens is a browser tool for multiplex serology. For each antigen it fits a
two-component Bayesian mixture model (a seronegative group and a seropositive group),
estimates seroprevalence, checks that the sampler worked, and tests how much the
estimates depend on the priors.

It is a single HTML page: **everything runs in your browser and no data leave your
computer.** There is nothing to install and no server.

## What it does

| Step | What you do |
|---|---|
| 1. Data | Load SerolyzeR output (a `sample_name` column and one column per antigen) and, optionally, a metadata file. Filter samples and choose the log transform. |
| 2. Baseline priors | Inform the priors yourself (per antigen, shared, or from a prior table), or use priors suggested by the data: a peak and tail rule or a maximum-likelihood mixture fit (EM). Set whether each antigen's prevalence is expected to be below or above 50%. |
| 3. Fit checks | R-hat, bulk and tail ESS, divergences, and whether the chains agree, for every fit. |
| 4. Prevalence estimates | Prevalence with 95% credible intervals for each antigen, overall or by group. |
| 5. Prevalence by covariates | Let prevalence differ between cells defined by metadata (for example subdistrict × age group), with shared negative and positive distributions. Models are compared with ELPD-LOO. |
| 6. Sensitivity analysis | Refit each antigen under a grid of alternative priors, with the base model or any step 5 model, and see which priors move the estimate. |

Use **Save session** (top right) to download a copy of the page with your data,
settings and results inside it. Open that file later to pick up where you left off.

The sampler is a JavaScript implementation of the No-U-Turn sampler (as in Stan). The
equivalent Stan models are in [`stan/`](stan/) and can be downloaded from the page,
together with an R script to rerun any analysis in Stan.

## Input format

**Antibody data** (required): CSV or TSV, one row per sample.

```
sample_name,DENV1,DENV2,...,CHIKV
CHPI_1,1.07,0.67,...,3.68
```

An optional `plate_name` column is recognised. Every other numeric column is treated as
an antigen.

**Metadata** (optional): CSV or TSV with one column of sample names matching the data
file (any name; the page finds it) and any number of covariate columns.

```
sample_name,upazila,ageG
CHPI_1,Chapai Nawabganj Sadar Upazila,5
```

See [`examples/`](examples/) for complete files. The page also opens with this example
loaded.

## Publish it on GitHub Pages

The repository is [github.com/Leacavalli/SeroPrevaLens](https://github.com/Leacavalli/SeroPrevaLens).

1. Push the files to the `main` branch:
   ```bash
   git add .
   git commit -m "Add SeroPrevaLens"
   git push
   ```
2. On GitHub, open **Settings → Pages**. Under **Build and deployment**, set **Source**
   to *Deploy from a branch*, choose **main** and **/ (root)**, and save.
3. After a minute or two the tool is live at
   <https://leacavalli.github.io/SeroPrevaLens/>. On a free account the repository
   must be public for Pages to work.

To update the tool, replace `index.html` and push again.

## Use it offline

Download `index.html` and open it in a recent Chrome, Edge, Firefox or Safari. The logo
fonts load from Google Fonts when online; without a connection the page falls back to
system fonts and everything else works.

## Repository contents

```
index.html                  the tool (self-contained)
assets/logo.svg             logo, light background
assets/logo-dark.svg        logo, dark background
assets/favicon.svg          icon without labels
stan/Mixture1D_priors.stan  1D mixture with priors and limits passed as data
stan/Mixture1D_cells.stan   one prevalence per cell, shared components (step 5)
examples/                   example antibody data and metadata (see below)
CITATION.cff                how to cite this tool
LICENSE                     licence for the code
.nojekyll                   tells GitHub Pages to serve the files as they are
```

## Example data and default priors

The example data and the default (Table 1) priors come from:

O'Driscoll M, Hozé N, Lefrancq N, Ribeiro dos Santos G, Hoinard D, Rahman MZ, Paul KK,
Titu AMN, Alam MS, Hossain ME, Vanhomwegen J, Cauchemez S, Gurley ES, Salje H.
Epidemiological and antigenic inferences from serological cross-reactivity among
arboviruses. *Science Translational Medicine* 17, eads8680 (2025).
[doi:10.1126/scitranslmed.ads8680](https://doi.org/10.1126/scitranslmed.ads8680).
Data and code: [doi:10.5281/zenodo.17340210](https://doi.org/10.5281/zenodo.17340210)
and [github.com/meganodris/MultiSero](https://github.com/meganodris/MultiSero).

These are anonymised measurements from a 2014 serosurvey of 1,453 people in Chapai
Nawabganj, Bangladesh. They are redistributed here under the terms of the original
deposit. The licence in [`LICENSE`](LICENSE) covers the code, not these data.

With the Table 1 priors the page reproduces the paper's reported 1D estimates (for
example DENV3 10.6%, 95% CrI 8.9–12.4%). The page itself, its sampler, the sensitivity
analysis and the data-suggested priors are an independent implementation, not the
authors' code.

## Citation

If you use SeroPrevaLens, please cite it (see [`CITATION.cff`](CITATION.cff); GitHub
shows a "Cite this repository" button) and the O'Driscoll et al. (2025) paper when you
use the Table 1 priors or the 1D model as described there.

## Licence

Code: MIT, see [`LICENSE`](LICENSE).
