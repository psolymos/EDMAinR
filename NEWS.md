# Version 0.3-0, 2023-08-21

* `edma_scale` function contributed by Kevin M. Middleton (PR #11).

# Version 0.2-1, 2021-12-07

* `read_xyz` don't change landmark names to be syntactically valid.

# Version 0.2-0, 2021-11-04

* Separate `edma_fdm_report` and `edma_gdm_report` functions.
* Renaming `Tobs_*` functions to `global_*` to avoid confusion with T-test.
* Updated documentation.

# Version 0.1-11, 2021-08-17

* Added `edma_report()` to reproduce WinEDMA output.
* Fixed CI calculation: `ref_denom` argument should not flip
  numerator and denominator.
* 2D and 3D plot revised to color edges and label landmarks.

# Version 0.1-10, 2021-07-19

* Bug fix: FDM methods had issues with B=0.
* Added landmark names to 2D plots, use `labels=TRUE`.
* 3D plots now preserve XYZ aspect ration of 1.
* `signif_only=FALSE` allows plotting the top/bottom percentage 
  of linear distances (edges) in 2D/3D plots.

# Version 0.1-8, 2020-12-19

* Added new function `gpa_fit` to estimate mean form based on GPA.

# Version 0.1-7, 2020-11-27

* Renaming `T_*` functions to `Tobs_*` to avoid confusion with T-test.
* Observed statistic not included into the null distribution for global testing.
* Added interactive plotly graphics to vignettes.

# Version 0.1-6, 2020-11-20

* Shape difference matrix (SDM) calculation added: `edma_sdm`.
* `edma_sdm` can be set to assume same size, use TLS, or a particular edge.
* CI calculation (local FDM, GM, GDM, SDM testing) calculations use
  2-sample bootstrap.
* Global FDM, GM, GDM testing (T-test) uses the mixed or the reference
  bootstrap, global Z-test for FDM is based on 2-sample bootstrap.

# Version 0.1-4, 2020-08-23

* `.edma_fit_np` is more preformant (original implementation
  retained as `.edma_fit_np_old` for comparison);
  it also exposes mean/variance for distances for future uses.
* `edma_simulate_data` is also much faster now.
* Nonparametric EDMA fid and FDM analysis can use multiple cores.
* Added `write.xyz` function to write EDMA data into xyz format.
* Added 4 data sets for growth difference analysis.

# Version 0.1-3, 2020-07-05

* Small fixes to parametric test helper functions.
* Miscellaneous function to print/visualize pattern matrices.
* Exposed `combine_data` and `combine_data4` functions.
* Clustering method now can be changed by the user.
* `as.edma_data` method added to turn a 3D array (a common
  morphometrics data format) to EDMA data objects.
* Added `edma_colors` and `plot_edma_colors` for
  manipulating color palettes set via `getOption("edma_options")`.
* Extensive updates to the docs and vignettes.

# Version 0.1-2, 2020-04-30

* Parametric estimation tested.

# Version 0.1-1, 2020-02-16

* Global test p-values updated.
* Package website made with {pkgdown}.

# Version 0.1-0, 2019-12-27

* EDMA package command line interface finalized.
