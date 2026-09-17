# Methods

## Public source

NCBI Gene Expression Omnibus accession **GSE58331**.

GEO describes:
- organism: *Homo sapiens*
- platform: GPL570, Affymetrix Human Genome U133 Plus 2.0
- processed sample `VALUE`: log2 RMA signal
- tissues: anterior orbit and lacrimal gland
- diagnoses: TED, normal controls and several other inflammatory orbital diseases
- some subjects with multiple samples / technical replicates

## Primary cohort

Only titles beginning:
- `TED Anterior Orbit`
- `Normal Anterior Orbit`

are included in the primary model.

Biological subject is parsed from disease class plus the first integer following the tissue label. This deliberately groups repeated samples from the same subject together before validation.

## Preprocessing

1. GEO submitter-processed log2 RMA `VALUE`
2. GPL570 probe → gene symbol annotation
3. median across probes mapped to the same gene symbol
4. mean across repeated / multi-biopsy samples assigned to the same biological subject
5. zero-variance genes removed before model fitting

No outcome-based gene filtering is performed outside cross-validation.

## QC and outliers

The primary model retains all eligible biological subjects.

Two descriptive QC signals are calculated:
- robust Mahalanobis distance in the first up-to-five PCA components
- each subject's median correlation to all other subjects

A subject is flagged if:
- PCA distance exceeds the 99th-percentile chi-square reference, or
- median-correlation robust z-score is < -3.5

These are **review flags**, not automatic exclusion rules. A secondary sensitivity analysis repeats model validation with flagged subjects removed.

For subjects with multiple source samples, pairwise repeat correlations are also saved.

## Model

Within each outer training fold:

1. median imputation
2. `SelectKBest(f_classif)`
3. standardization
4. class-balanced L2 logistic regression

Inner CV tunes:
- number of selected genes: 5 / 10 / 20 / 40 / 80
- logistic regularization C: 0.01 / 0.1 / 1 / 10

Primary validation:
- repeated stratified 5-fold outer CV
- 10 repeats
- 3-fold inner tuning

Every biological subject receives one held-out prediction per repeat. The dashboard's subject score is the mean of those repeated held-out predictions.

## Performance outputs

- AUROC
- average precision
- Brier score
- bootstrap 95% interval around subject-level OOF AUROC
- descriptive calibration bins
- threshold-specific sensitivity, specificity, PPV, NPV and balanced accuracy
- feature-selection frequency
- subject score SD across repeated outer CV

## Leakage stress test

A deliberately naive sample-level validation is also run. It does not group repeated samples from the same subject and is **not** reported as the valid primary estimate.

Its only purpose is to demonstrate how repeated biological material can inflate apparent performance when allowed to cross the train/test boundary.

## Tissue-domain stress test

A final anterior-orbit model is applied to TED / normal lacrimal-gland subjects from the same GEO series.

This evaluates tissue transportability, not independent clinical validation. If patient identifiers overlap across anatomical sites, the result should be interpreted as a partially paired cross-tissue projection rather than external validation.
