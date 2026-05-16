# MRBEE CLI

Command-line interface for running [MRBEE](../README.md) multivariable Mendelian Randomization analyses. Analysis parameters are supplied through a YAML config file; no R scripting required.

## Requirements

[pixi](https://pixi.sh) must be installed.

## Setup

```sh
# 1. install conda-forge / bioconda packages (R, plink)
pixi install

# 2. install CRAN-only dependencies (FDRestimation) and the MRBEE package itself
pixi run setup
```

`setup` only needs to be run once (or again after updating `../`).

## Config file

Copy `analyses/template.yaml` and fill in the paths and column names for your GWAS files.

```yaml
ld: data/1000G/1kg_phase3_EUR_only   # PLINK bfile prefix — no .bed/.bim/.fam suffix

outcome:
  file: path/to/outcome.tsv
  rsid_column: SNP
  effect_allele_column: A1
  non_effect_allele_column: A2
  effect_size_column: BETA
  standard_error_column: SE

exposures:
  ldl:
    file: path/to/ldl.tsv
    rsid_column: SNP
    effect_allele_column: A1
    non_effect_allele_column: A2
    effect_size_column: BETA
    standard_error_column: SE
  hdl:
    file: path/to/hdl.tsv
    rsid_column: rsid
    effect_allele_column: effect_allele
    non_effect_allele_column: other_allele
    effect_size_column: b
    standard_error_column: se
```

Each exposure key (e.g. `ldl`, `hdl`) becomes the exposure label in the output. Column names do not need to match across files. GWAS files can be in any delimiter-separated format readable by `data.table::fread` (tab, comma, space, gzipped, etc.).

## Usage

```sh
pixi run mrbee -- -c analyses/my_analysis.yaml
```

With options:

```sh
pixi run mrbee -- \
  -c analyses/my_analysis.yaml \
  -o results/my_analysis.csv \
  --var-est ordinal \
  --pv-thres 0.01
```

| Flag | Default | Description |
|---|---|---|
| `-c`, `--config` | *(required)* | Path to the analysis config YAML |
| `-o`, `--out` | `mrbee_results.csv` | Output CSV path |
| `--var-est` | `variance` | Residual variance estimator: `variance`, `robust`, or `ordinal` |
| `--pv-thres` | `0.05` | P-value threshold for the pleiotropy test |
| `--no-fdr` | off | Disable FDR correction for the pleiotropy test |

## Pipeline

```
[1/6] Read config
[2/6] Load GWAS files
[3/6] Harmonize alleles (HapMap3 reference panel)
[4/6] Estimate error covariance matrix (GWAS-insignificant SNPs)
[5/6] Select IVs — joint test (MRBEE::Joint.test) + PLINK C+T (p1 = 5e-8, r2 = 0.01, kb = 500)
[6/6] Estimate causal effects (MRBEE::MRBEE.IMRP)
```

The error covariance matrix (`Rxy`) is estimated from all harmonized SNPs before IV selection, so insignificant SNPs contribute to the covariance estimate even if they are not used as IVs. IV selection runs `MRBEE::Joint.test` on the exposure Z-scores to get a joint p-value per SNP, then uses PLINK to LD-clump those p-values and retain independent genome-wide significant SNPs.

## Output

A CSV with one row per exposure:

| Column | Description |
|---|---|
| `exposure` | Exposure label from the config |
| `estimate` | Causal effect estimate (beta) |
| `std_error` | Standard error of the estimate |
| `z_stat` | Wald z-statistic |
| `p_value` | Two-sided p-value |
