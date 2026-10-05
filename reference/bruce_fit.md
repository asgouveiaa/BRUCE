# Fit the model for one time window

Fit the model for one time window

## Usage

``` r
bruce_fit(
  d,
  mesh,
  spec = bruce_spec(),
  window = NULL,
  theta_init = NULL,
  fix_theta = FALSE,
  keep_config = FALSE,
  cpo = FALSE,
  drop_stations = NULL,
  verbose = FALSE
)
```

## Arguments

- d:

  A `bruce_data` object (usually after
  [`bruce_covariates()`](https://asgouveiaa.github.io/bruce/reference/bruce_covariates.md)).

- mesh:

  Mesh from
  [`bruce_mesh()`](https://asgouveiaa.github.io/bruce/reference/bruce_mesh.md).

- spec:

  Model specification from
  [`bruce_spec()`](https://asgouveiaa.github.io/bruce/reference/bruce_spec.md).

- window:

  Two dates (first and last day). `NULL` uses the full range of `d$obs`.
  Windows are fitted independently; one month is the tested scale.

- theta_init:

  Optional starting mode for the hyperparameters (on INLA's internal
  scale, e.g. `other_fit$inla$mode$theta` from a stable fit with the
  same structure).

- fix_theta:

  Keep the hyperparameters fixed at `theta_init` instead of estimating
  them: INLA then computes a single Gaussian approximation of the latent
  field, which is much faster. Meant for operational maps with
  hyperparameters learned from a previous fit, and for fast
  cross-validation. Hyperparameter uncertainty is not propagated, and
  `theta_init` must come from a model with the same structure (same
  networks with estimated noise, same optional terms).

- keep_config:

  Keep the latent configuration so that posterior samples can be drawn
  (needed by `bruce_predict(n_samples > 0)`).

- cpo:

  Compute CPO/PIT.

- drop_stations:

  Station ids excluded from fitting (used by
  [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md)).

- verbose:

  Passed to INLA.

## Value

An object of class `bruce_fit`.

## Details

The linear predictor for an observation of network `k` at station `i`,
location `s` and day `t` is \$\$\eta = \beta_0 + \delta_k +
x(s)^\top\beta + u_i + d_t + w(s, t)\$\$ where \\\delta_k\\ is the
network offset (only for `bias = TRUE`), \\u_i\\ the iid station effect,
\\d_t\\ the city-wide daily AR(1) term and \\w(s,t)\\ a Matern SPDE
field with AR(1) dependence between days. Each network, and the
background, has its own Gaussian likelihood. Background cells share
\\\beta_0\\, covariates, \\d_t\\ and \\w\\ but no offset or station
effect.

After fitting, the field is checked for collapse (see
[`bruce_spec()`](https://asgouveiaa.github.io/bruce/reference/bruce_spec.md));
every attempt is recorded in `fit$attempts` and the outcome in
`fit$status`.

The mesh resolution is also checked against the estimated range: the
typical edge of the triangles that contain the data is compared with the
posterior mean range (`fit$mesh_resolution`). When the range is less
than five times that edge, a warning is issued: the field is poorly
resolved and the range estimate can change with the triangulation (in
the Rio de Janeiro application, two meshes of about 1,000 vertices built
with the same settings by different fmesher versions gave ranges of 7.7
and 4.1 km).
