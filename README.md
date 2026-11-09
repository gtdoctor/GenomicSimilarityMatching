# PCA Genomic Similarity Matching

GOAL: match controls to cases based on distances in the space defined by principal components of a standardised genotype matrix of all cases and controls. Remove cases that cannot be adequately matched. Output independent strata — each retained case with its ncmin–nc matched controls, no control shared between strata.

## INPUTS:

Requires:
**1. PCA eigenvalues** (`PCA.EIGENVAL`) one value per line, no header.
Must be eigenvalues (variance per PC), not singular values. PLINK's `.eigenval` is correct as-is; `bigsnpr`'s `$d` holds singular values and must be squared first.
**2. PC scores** (`PCA.EIGENVEC`) with columns `FID`, `IID`, `PC1`…`PC10`. Headerless, or a header line beginning `#`.
Must be either:

- principal component scores - samples projected onto PC axes. This is what `bigsnpr` returns via `predict()`, the R projection function; or
- eigenvalues (e.g. the output of plink/plink2 --pca - these are unit vector, standard-deviation normalised (i.e. whitened principal component scores)
   **3. Sample file** (`SAMPLEFILE`) one row per sample, columns addressed by index:

* `SAMPLEFILE_IDCOL` (default 2) — sample ID, matching `IID` in the scores file
* `SAMPLEFILE_PHENOCOL` (default 6) — `1` = control, `2` = case
* `SAMPLEFILE_BATCHCOL` (optional) — batch/chip label, to restrict which controls each case may match

## OUTPUTS:

__Part A__ — `<stem>_<n>pcs[_scaled]_<metric>.Rdata`, holding the matrix `ccdistances` (cases × controls).

__Part B__ — a directory `<stem>_maxt<>mint<>_minm<>maxm<>_<date>/` containing three headerless .tsvs:

* `filteredmatches.*` — the strata: `control_id`, `case_id`, `distance`; ordered by case then distance. One row per case–control pair; each control appears exactly once across the whole file.
* `kept.*` — IIDs of retained cases and controls, preceded by `#` comment lines recording the parameters and summary counts.
* `removed.*` — IIDs of dropped cases and controls.

**Part C** — three JPEGs of PC1 against PC2–PC5:

* `plot.allsamples.*` — all cases vs all controls, before matching
* `plot.splitsamples.*` — cases and controls split by kept/removed
* `plot.onlykept.*` — retained cases and controls only

## PART A: Distance Matrix

Creates a matrix of distances between each case and all controls.

### Parameters:

* `DATADIR`  = directory of input data
* `OUTDIR` =  directory of output
* `PCA.EIGENVAL` as above
* `PCA.EIGENVEC` as above
* `SAMPLEFILE_IDCOL` as above
* `SAMPLEFILE_PHENOCOL` as above
* `SAMPLEFILE` = file containing IDs and phenotype code
* `FILEOUTSTEM` = user-defined part  of output names
* `SMALLSAMPLE` = FALSE; TRUE subsets to the first 500 cases for testing.
* `PCmax` = Number of principal components to include (max 10)
* `euclidean` = TRUE for Euclidean distance, FALSE for Manhattan (sum of absolute differences along each PC axis).
* `eigvalweight` = FALSE (default). TRUE weights each PC score/eigenvector by its eigenvalue as a proportion of the total, w_k = λ_k / Σλ.
* `scale`  = normalise each PC score/eigenvector by its standard-deviation (i.e. whiten). This converts bigsnpr-output PC scores to plink --pca output eigenvectors.
* `nchunks` = 10 (default) number of chunks for parallelisation of large datasets

For PLINK --pca style eigenvectors, `scale` should be set to FALSE. (In fact it has no effect since every PC already has the same standard deviation.)  `eigvalweight = TRUE` converts such eigenvetors to principle component scores, equivalent to the bigsnpr output. Conversely, bigsnpr output with  `scale=TRUE` and `eigvalweight` =TRUE is the same as if both those parameters = FALSE.

## PART B: Matching

* Based on a distance matrix between cases and controls, finds between  `ncmin` and `nc` unique controls for each case within a maximum distance.
* Cases are removed if fewer than `ncmin` controls can be found, along with  any controls provisionally assigned to them.
* Each control is used at most once, so the retained cases and their controls  form independent strata.
* Cases are removed if they match to fewer than "ncmin" controls.

### Parameters:

### Parameters:

* `wdt` = distance matrix from PART A
* `datadir` / `outdirp` = input directory / parent of the output directory
* `SAMPLEFILE`, `SAMPLEFILE_IDCOL`, `SAMPLEFILE_PHENOCOL` as in PART A
* `fileoutstem` = user-defined part of output names
* `nc` = maximum controls per case
* `ncmin` = minimum acceptable controls per case; cases below this are excluded
* `maxthresh` = upper distance cutoff, as a multiple of the median   admissible case–control distance (1 = the median itself; smaller is tighter)
* `minthresh` = lower distance cutoff, as a quantile of the admissible   distance distribution (0 = no floor). Note this is a probability, not a multiple of the median.  A floor can be used to avoid matching cryptically related pairs. (To Do: Thresholding to be aligned.)
* `SAMPLEFILE_BATCHCOL` = optional column holding a batch/chip label,  restricting which controls each case may match
* `MATCHSELF` = TRUE allows a case to match controls from its own batch
* `COMMONBATCH` = batch (or batches) of controls matchable by every case regardless of their own batch. With `MATCHSELF=FALSE`, only these controls are available and `COMMONBATCH` is required.
* `FORBIDDEN` = optional file of case–control pairs to exclude
* `useMILP` = TRUE adds stage 2 below. Requires R >= 4.5.1.
* `STRATA_REASSIGN` = FALSE (default). TRUE adds stage 3 below.
* `REFINE_TOPK` = 20 (default). Candidate controls per case in the stage 3 graph (of only matchable cases and controls); only used when `STRATA_REASSIGN=TRUE`.

Both cutoffs are computed after batch and forbidden-pair exclusions, so they describe the distribution of permissible matches rather than all pairs.

### Further details of matching:

Three stages are applied, the last two optional.

1. **Localbestworst** (always)

   * Ignoring distances outside the cutoffs, each case's nearest available control is identified.
   * Among cases still short of `ncmin`, the case whose nearest control is furthest is assigned that control — worst-best first, so the hardest cases are served before the easy ones.
   * That control is withdrawn from the pool. Every case that had been pointing at it finds its next nearest. Repeat until no case can gain a control or all have `nc`.

2. **Min-cost-max-flow on deficit cases** (`useMILP=TRUE`)

   * Deficit cases are those with 0 < matches < `ncmin`. Taking those cases and the controls reserved to them, a binary MILP first maximises how many      reach `ncmin`; a min-cost LP then assigns controls to the cases it found.
   * Important for rescuing rarer-ancestry cases, where stage 1 can divide a small pool of controls between several cases without any reaching the
      minimum.

3. **Reassignment** (`STRATA_REASSIGN=TRUE`)

   * This stage reallocates the already-selected controls among the already-retained cases to minimise average distances, each control used exactly
      once and each case holding between `ncmin` and `nc`.
   * It introduces no new controls and drops no cases: the kept sets are unchanged and only the pairings move.
   * Only relevant for strata-aware analysis.

## COMMENTS

* Results are sensitive to the method chosen

## PART C: Plotting

Run from the command line with three positional arguments:

    Rscript part_C.R <PC scores> <sample file> <kept file>

* `<PC scores>` — the same scores/eigenvector file used in PART A
* `<sample file>` — the same sample file used in PARTS A and B
* `<kept file>` — the `kept.*.tsv` written by PART B. Use the full path.

### Parameters

* `PCMAX` = 10. Must equal the number of PC columns in the scores file.
* `NORMPLOT` = "y". Scales each PC to unit length for plotting only, so all PCs share an axis scale. Set "n" to plot raw values, where PC1 dominates and the later PCs are compressed to a thin band.

## Developers

Ancestry matching algorithm first developedby AP Levine. https://github.com/APLevine
Updated 2024-6 by GT Doctor UCL Centre for Genetics and Genomics https://github.com/gtdoctor