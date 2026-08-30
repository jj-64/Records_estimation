# records: Record-count statistics

This package provides utilities to compute expected number of records, variances,
and exact distributions of the number of records for several stochastic models:

- IID (independent identical)
- DTRW (discrete-time random walk)
- LDM (linear drift model)
- YNM (Yule–Simon / preferential attachment family)

The package uses a registry-based dispatcher. See `rec_count_stats()` for unified access.
