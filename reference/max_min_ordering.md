# Maximum minimum distance ordering

Returns the indices of an exact maximum-minimum distance ordering. This
is similar to the
[`order_maxmin_exact`](https://rdrr.io/pkg/GPvecchia/man/order_maxmin_exact.html)
function. The main difference is that it allows the user to specify a
distance function.

## Usage

``` r
max_min_ordering(locs, dist.func = NULL)
```

## Arguments

- locs:

  A matrix with one row per location and any number of columns (x, y,
  time, etc.).

- dist.func:

  Any distance function with a signature of dist(query_location,
  locations_matrix). Defaults to Euclidean distance.

## Value

A vector of indices giving the ordering. Element *i* of this vector is
the index of the *i*th location.

## Details

For Euclidean distance, this function will return the same results as
[`order_maxmin_exact`](https://rdrr.io/pkg/GPvecchia/man/order_maxmin_exact.html),
but will be much slower for large data sets.
[`order_maxmin_exact`](https://rdrr.io/pkg/GPvecchia/man/order_maxmin_exact.html)
should be used instead of this function when the distance function is
Euclidean.

## References

- Katzfuss, M., and Guinness, J. "A general framework for Vecchia
  approximations of Gaussian processes", Statistical Science (2021)
  36(1):124-141.

- Guiness, J. "Permutation methods for sharpening Gaussian process
  approximations", Technometrics (2018) 60(4):415-429.

## See also

[`order_maxmin_exact`](https://rdrr.io/pkg/GPvecchia/man/order_maxmin_exact.html)

## Examples

``` r
data(weather)
locs <- weather[,3:4]
max_min_ordering(locs)
#>   [1] 131  97   1  45  82  69 109 155  39   5  81  27 129  89  99  14  29  32
#>  [19]  88   4  52 148 125  58 111 124 113  13  48  44 130  92  84  51  79 100
#>  [37]  41  86  34 116 122  16  11  40 112  57  60 103  26  98 133 108  49   2
#>  [55] 134  95  20 128  65 120  75  54 138 142 153 105  21  56 137  10  61  93
#>  [73] 152  83  47 156   3  24 119  68 140   9  50 135  96 146  91  42 149  37
#>  [91] 110 106  19 141  62 123 151  66  80  38 127 145  53 118  23 104 136  63
#> [109] 114  94  59  55  25 126  87  30 157  18  70 101  64 139 115 107  33   7
#> [127]  71  90 144 117  17  31 121 143  78  73  22 154 102  43  36 132  67  85
#> [145]  77  35 147  74 150  46  72  15   8  76  12  28   6
```
