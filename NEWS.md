# v 0.2.5

* Added two functions, `information_volume` and `correspondence` to calculate
weights based on the amount of things that are non-zero in both things when
doing pairwise correlation.

* Added logical argument `weight` to `pairwise_correlation` to weight the correlations. If `weight = TRUE`, the diagonal will not be **1** anymore, but instead will reflect how many features out of the total are in that sample.

# v 0.2.3

* A bug was discovered in `median_correlations` that meant the wrong sample ids
might be added to the output data, making detection of real problems difficult

# v 0.2.1

* `pairwise_correlation` now uses `cor` internally directly, whereas previously
it did a `for` loop to allow pairwise comparisons. This makes the correlations
3x faster.

* `count` has been removed from the list returned by `pairwise_correlation`

* new function `pairwise_correlation_count` to get the counts in each pairwise
comparison
 
# v 0.1.1

* Changed correlation function to return a list instead of a matrix. This
list contains the correlations (`cor`), counts in each correlation (`count`),
and which points passed the criteria (`keep`).
