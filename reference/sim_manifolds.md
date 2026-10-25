# Simulated manifolds

The Swiss roll and S-curve manifolds of Gashler, Ventura and Martinez
(2008), with their intrinsic (two-dimensional) coordinates, to test
dimensionality reduction methods.

## Usage

``` r
sim_swiss_roll(n = 2000, noise = 0)

sim_s_curve(n = 2000, noise = 0)
```

## Arguments

- n:

  The number of points.

- noise:

  The standard deviation of Gaussian noise added to the coordinates `x`,
  `y` and `z`.

## Value

A tibble with `n` rows and the columns `x`, `y`, `z` (the coordinates in
three dimensions) and `u`, `v` (the intrinsic coordinates). The points
are ordered by `u`.

## Details

For the Swiss roll, \\t_i = 8i/n + 2\\ for \\0 \le i \< n\\, and \$\$x_i
= t_i \sin(t_i), \quad y_i \sim U(-6, 6), \quad z_i = t_i \cos(t_i).\$\$
Its intrinsic coordinates are the arc length along the roll, \\u_i =
(\sinh^{-1}(t_i) + t_i \sqrt{t_i^2 + 1}) / 2\\, and \\v_i = y_i\\.

For the S-curve, \\t_i = (2.2i - 0.1)\pi / n\\, and \$\$x_i = t_i, \quad
y_i = \sin(t_i), \quad z_i \sim U(0, 2).\$\$ Its intrinsic coordinates
are the arc length along the curve, \\u_i = \int_0^{t_i}
\sqrt{\cos^2(w) + 1} \\ dw\\, and \\v_i = z_i\\. (The paper gives \\v_i
= y_i\\, a misprint since \\y_i\\ is a function of \\u_i\\.)

## References

Gashler, M., Ventura, D. and Martinez, T. (2008). Iterative non-linear
dimensionality reduction with manifold sculpting. In *Advances in Neural
Information Processing Systems 20*, 513–520.

## Examples

``` r
set.seed(1)
sim_swiss_roll(5)
#> # A tibble: 5 × 5
#>       x      y      z     u      v
#>   <dbl>  <dbl>  <dbl> <dbl>  <dbl>
#> 1  1.82 -2.81  -0.832  2.96 -2.81 
#> 2 -1.59 -1.53  -3.23   7.72 -1.53 
#> 3 -4.59  0.874  2.44  14.9   0.874
#> 4  3.36  4.90   5.91  24.7   4.90 
#> 5  7.18 -3.58  -4.36  36.9  -3.58 
sim_s_curve(5)
#> # A tibble: 5 × 5
#>         x       y     z       u     v
#>     <dbl>   <dbl> <dbl>   <dbl> <dbl>
#> 1 -0.0628 -0.0628 1.80  -0.0888 1.80 
#> 2  1.32    0.969  1.89   1.66   1.89 
#> 3  2.70    0.426  1.32   3.21   1.32 
#> 4  4.08   -0.809  1.26   5.07   1.26 
#> 5  5.47   -0.729  0.124  6.54   0.124
```
