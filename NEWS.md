# v 0.2.1

* `pairwise_correlation` now uses `cor` internally directly, whereas previously
it did a `for` loop to allow pairwise comparisons. This makes the correlations
3x faster.

* `count` has been removed from the list returned by `pairwise_correlation`

* new function `pairwise_correlation_count` to get the counts in each pairwise
comparison

* new function
 
# v 0.1.1

* Changed correlation function to return a list instead of a matrix. This
list contains the correlations (`cor`), counts in each correlation (`count`),
and which points passed the criteria (`keep`).
