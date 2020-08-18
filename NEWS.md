# Version 0.1-4, 2020-08-17

* `.edma_fit_np` is more preformant (original implementation
  retained as `.edma_fit_np_old` for comparison);
  it also exposes mean/variance for distances for future uses.
* `edma_simulate_data` is also much faster now.

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
