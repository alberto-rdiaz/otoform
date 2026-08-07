
# otoform

Extraction and normalisation of otolith contours from photographs, for
morphometric shape analysis via à trous wavelet decomposition of
centroid-distance signatures.

**Status:** early development. The API may change without notice.

## Installation

```r
# install.packages("remotes")
remotes::install_github("alberto-rdiaz/otoform")
```

## Image requirements

The pipeline assumes the following. Images that do not meet these
conditions will produce incorrect results **without raising an error**.

- **Dark, uniform background**, as free of dust as possible.
- **A single otolith**, centred in the frame.
- **Rostrum pointing right.** The starting point of the contour is defined
  as the rightmost point after rotation onto the principal axis.
- **Consistent resolution** across the collection.

## Usage

```r
library(otoform)

res <- process_directory("path/to/images")
write_results(res, "path/to/output")
```

This writes three files:

| File | Contents |
|---|---|
| `otoform_results.rds` | Complete results, including raw contours |
| `otoform_coords.txt` | Normalised coordinates, long format |
| `otoform_log.txt` | One row per image: status, parameters, diagnostics |

Always check the log before analysing the output. Images that fail are
recorded there rather than interrupting the run.

## Pipeline

1. Greyscale conversion and binarisation
2. Contour extraction by Moore boundary tracing, forced counter-clockwise
3. Resampling to a uniform arc-length grid
4. Smoothing with a binomial kernel
5. Centring and rotation onto the principal axis
6. Sub-pixel rostrum anchoring and resampling to a fixed number of points

The final contour has the rostrum at index 1, to sub-pixel accuracy. This
makes contours from different specimens directly comparable point by point.

## License

MIT
