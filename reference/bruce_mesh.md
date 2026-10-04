# Build the SPDE mesh

Triangulates around all observation locations (stations and background
cells) with a non-convex hull boundary.

## Usage

``` r
bruce_mesh(d, max_edge = NULL, cutoff = NULL, offset = NULL, convex = -0.1)
```

## Arguments

- d:

  A `bruce_data` object.

- max_edge:

  Length-2 vector: maximum triangle edge inside and outside the hull
  (metres).

- cutoff:

  Minimum distance between mesh vertices (metres).

- offset:

  Length-2 vector: inner and outer extension (metres).

- convex:

  Convexity parameter passed to
  [`fmesher::fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.html).

## Value

An `fm_mesh_2d` object with attribute `bruce_nn` (the station
nearest-neighbour median, used by the default range prior).

## Details

When left `NULL`, sizes are derived from the data: the median
nearest-neighbour distance between stations (`nn`) and the diagonal of
the bounding box (`L`). The inner edge is
`clamp(2 * nn, L / 40, L / 10)`, the outer edge 8 times that, the cutoff
0.6 times the inner edge, and the offsets `inner / 5` and `L / 2`. These
are heuristics; for a city-scale network they land in the same range as
the values used in the original Rio de Janeiro application
(`max_edge = c(5000, 40000)`, `cutoff = 3000`,
`offset = c(1000, 50000)`), but check the number of vertices (`mesh$n`)
before fitting a full window.
