# otoform (development version)

* `toBW()` removed. `binarize(method = "fixed", threshold = 127L)` produces
  a pixel-identical mask on the 10 reference images (0 differing pixels).
* `process_image()` and `process_directory()` gain `method` and `threshold`,
  passed on to `binarize()`. The default is now Otsu's method.
* Fixed the hint in `binarize()`'s range error: a cut at 0.5 corresponds to
  `threshold = 127`, not 128.
