# otoform 0.1.2

# otoform 0.1.1

# otoform 0.1.0

# otoform (development version)

* `toBW()` removed. `binarize(method = "fixed", threshold = 127L)` produces a
  pixel-identical mask on the 10 reference images (0 differing pixels).
* `process_image()` and `process_directory()` gain `method` and `threshold`,
  passed on to `binarize()`. The default is now Otsu's method.
* Fixed the hint in `binarize()`'s range error: a cut at 0.5 corresponds to
  `threshold = 127`, not 128.
* New `normalize_area()`: rescales a closed contour to unit area, keeping the
  original area as an attribute.
* New `centroid_distance()` and internal `polygon_moments()` /
  `polygon_centroid()`: centroid-distance signature referenced to the area
  centroid of the contour.
* `process_image()` also returns the signature `r`, the enclosed `area` in
  squared pixels, and `center`.
* New `process_contours()`: applies `atrous_b3()` to the signature of every
  element returned by `process_directory()`.
* New `keep_largest()` (via `mmand`): keeps the largest connected component of
  the mask. `Conte()` now starts at the first object pixel in raster order, so
  interior holes and background debris can no longer be traced.
* `process_image()` checks that the traced polygon covers its component
  (`min_area_ratio = 0.95`) and records `area_ratio` and `n_components`.
* `rotation()` centres on the area centroid and takes the principal axis from
  the polygon's second area moments, consistent with `centroid_distance()`;
  returns `center`, `anisotropy` and `status`.
* New `write_overlay()` and internal `to_image()`: PNG of the image with the
  contour, rostrum and centroid drawn; `process_image()` gains `overlay_file`.
* `process_directory()` writes `log.csv`, `results.rds` and optional overlays
  to `out_dir` (default `dir/otoform`), in blocks via `future.apply`; parallel
  when the user sets a future plan.
* `process_directory()` now runs both stages (`J` argument) and writes
  `contours.csv`, `signature.csv` and one `wavelet_*.csv` per level, plus
  `results.rds` and `wavelets.rds`. Export functions renamed `export_*`
  (`export_log`, `export_contours`, `export_signatures`, `export_wavelets`,
  `export_overlay`); `coords_long()`, `write_results()`, `build_log()` and
  `write_overlay()` removed.
* `J` is no longer a hard-coded default. New `max_levels(n)` applies the
  half-turn criterion (`2^(J+2) - 3 <= n/2`; 6 for n = 512, 7 for 1024);
  `atrous_b3()`, `process_contours()` and `process_directory()` take
  `J = NULL` to use it, and explicit values above the full-contour bound are
  rejected. `atrous_b3()` now rejects coordinate matrices and is vectorised:
  identical output, about 8x faster.
  * New `measure_contour()` and `export_measurements()`: area, perimeter,
  length and width along the principal axes, roundness, aspect and
  rectangularity from the aligned contour, in pixels or, with `cal` (px/mm),
  in millimetres. `process_directory()` writes `measurements.csv`.
* Sampling decisions fixed and documented: `n_points = 512` (set by the
  lowest-resolution images; finer levels are below the pixel-noise floor),
  `k = 32` (anti-aliasing at half a sampling step), `J = max_levels(N)`.
* process_image() records width and height of the processed image in the log. 
  Resolution sensitivity test documented in the vignette.
