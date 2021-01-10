# MSbox2

Common mass spectrometry tools for minimalist MS users.

# Functions

1. Check element isotopes

examples:

```r
E_iso('C') # element symbol, case insensitive
E_iso('Carbon') # element full name, case insensitive
E_iso('carBon') # element full name, case insensitive
```

2. Get extra molecular information based on compound name

If the queried information is not found for a compounds, it will assign "unknown" to that compound information:
```r
describe('malic acid') # get formula by default
describe(c('malic acid', 'citric acid', 'tartaric acid'), representation = "smiles") # get smiles
```

3. Calculate monoisitopic mass

example:

```r
mass('c7h6O1') # case insensitive: elements are seperated by numbers.
mass('C7H6O', caseSensitive = T) # case sensitive: elements are seperated by upper case letters. The number of the element can be missing if it is 1.
mass(c('K1', 'C5H8', 'nA1')) # vector input
mass(c('K1', 'C5H8', 'Na'), caseSensitive = T) # vector input
```

4. Calculate exact m/z values

example:

```r
mz('c7h6O1', z = -1) # case insensitive: elements are seperated by numbers.
mz('C7H6O', z = -1, caseSensitive = T) # case sensitive: elements are seperated by upper case letters. The number of the element can be missing if it is 1.
mz(c('C5H8', 'c3h3o1'), z = 1) # vector input
mz(c('C5H8', 'C3H3O'), z = 1, caseSensitive = T) # vector input
```

5. Calculate the mass accuracy of measured m/z

examples:

```r
ppm(155.03383, 155.03388) # with m/z value
ppm(155.03383, .03388) # lazy input when the integer parts of m and t are the same
ppm(155.03383, .03388, lazy = F) # lazy input disabled
ppm(155.03384, mz('C7H7O4', z = 1)) # with ion formula
```

6. Calculate isotope labelled molecular mass

example

```r
Iso_mass(F = 'C7H6O4', iso = '[13]C2[2]H3') # Two 13C and three 2H are labled. Case insensitive.
```

7. Calculate isotope labelled m/z

example

```r
Iso_mz(F = 'C7H6O4', iso = '[13]C2[2]H3', z = 1) # Two 13C and three 2H are labled. Case insensitive.
```

8. Check if an m/z value originates from possible contaminant

examples

```r
contam(33.0335, ppm = 10, mode = '+')
contam(44.998, ppm = 10, mode = '-')
```

9. Calculate common adduct ions in pos or neg ion mode

examples

```r
adduct('C1H4', mode = '-') # case insensitive
adduct('C1H4', mode = '+') # case insensitive
```

10. Annotate the m/z values according to m/z value matching

```r
what(1034.556, mode = "+", ppm = 3) # single m/z value

mzs <- rep(133.014, 300)
what(mzs, "-") # multiple m/z values, default ppm = 5
```
