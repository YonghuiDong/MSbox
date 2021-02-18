# MSbox

[![CRAN status](http://www.r-pkg.org/badges/version/MSbox)](https://cran.r-project.org/package=MSbox) 
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/MSbox)](https://cran.r-project.org/package=MSbox)

Common mass spectrometry and statistical tools.

# Installation

```r
install.packages('MSbox')
```

# Usage

## 1. Mass Spectrometry Widgets
-----------------------

### 1. Check element isotopes

examples:

```r
E_iso('C') # element symbol, case insensitive
E_iso('Carbon') # element full name, case insensitive
E_iso('carBon') # element full name, case insensitive
```

### 2. Get extra molecular information based on compound name

If the queried information is not found for a compounds, it will assign "unknown" to that compound information:
```r
describe('malic acid') # get formula by default
describe(c('malic acid', 'citric acid', 'tartaric acid'), representation = "smiles") # get smiles
```

### 3. Calculate monoisitopic mass

It accepts two types of inputs: 

(1) **Standards elemental composition**, for instance, `C7H6O`, `C7H6Na`. Here each element is distinguished by Capital letters, i.e., sodium should be written as Na, not NA, na or nA here. Since there is only one sodium in the formula `C7H7Na`, you don't have to write 1 after `Na`. 

(2) **User friendly elemental composition**, for instance,`c7H6O1`, `C7H6NA1`. Here each element is distinguished by the number of the element, therefore, i.e, sodium can be written as `Na`, `NA` or `na`. But the number the sodium element should be clearlly stated in the formula even if there is only one sodium. 

example for **Standards elemental composition**:

```r
mass("C7H6O", caseSensitive = T)
mass(c("C7H6O", "C7H6Na"), caseSensitive = T) # vector input
```

example for **User friendly elemental composition**:

```r
mass("c7h6O1")
mass(c("c7h6o1", "C7H6NA1")) # vector input
```

### 4. Calculate theoretical m/z value

It accepts two types of inputs: 

(1) **Standards elemental composition**, for instance, `C7H6O`, `C7H6Na`. Here each element is distinguished by Capital letters, i.e., sodium should be written as Na, not NA, na or nA here. Since there is only one sodium in the formula `C7H7Na`, you don't have to write 1 after `Na`. 

(2) **User friendly elemental composition**, for instance,`c7H6O1`, `C7H6NA1`. Here each element is distinguished by the number of the element, therefore, i.e, sodium can be written as `Na`, `NA` or `na`. But the number the sodium element should be clearlly stated in the formula even if there is only one sodium. 

example for **Standards elemental composition**

```r
mz("C7H7O", 1, caseSensitive = T) # [M+H]+, positive ion mode, charge z = 1
mz("C7H5O", -1, caseSensitive = T) # [M-H]-, negative ion mode, charge z = 1
mz(c("C7H6O", "C7H6Na"), 1, caseSensitive = T) # vector input
```

example for **User friendly elemental composition**:

```r
mz("c7h7O1") # [M+H]+, positive ion mode, charge z = 1
mz("c7H5o", -1) # [M-H]-, negative ion mode, charge z = 1
mz(c("c7h6o1", "C7H6NA1"), 1) # vector input
```

### 5. Calculate the mass accuracy of measured m/z ppm

example:

```r
ppm(155.03383, 155.03388) # standard way
ppm(155.03383, .03388) # lazy input when the integer parts of m and t are the same
ppm(155.03383, .03388, lazy = F) # lazy input disabled
ppm(155.03384, mz('C7H7O4', z = 1)) # with ion formula
```


### 6. Calculate isotope labelled molecular mass

example

```r
Iso_mass(F = 'C7H6O4', iso = '[13]C2[2]H3') # Two 13C and three 2H are labled. Case insensitive.
```

### 7. Calculate isotope labelled m/z

example

```r
Iso_mz(F = 'C7H6O4', iso = '[13]C2[2]H3', z = 1) # Two 13C and three 2H are labled. Case insensitive.
```

### 8. Check if an m/z value originates from possible contaminant

examples

```r
contam(33.0335, ppm = 10, mode = '+')
contam(44.998, ppm = 10, mode = '-')
```

### 9. Calculate common adduct ions in pos or neg ion mode

examples

```r
adduct('C1H4', mode = '-') # case insensitive
adduct('C1H4', mode = '+') # case insensitive
```

### 10. Search m/z in HMDB or KEGG database

```r
what(1034.556, mode = "+", ppm = 3) # single m/z value in HMDB database (default)
what(1034.556, mode = "+", ppm = 3, useDB = "KEGG") # single m/z value in KEGG database
what(c(133.014, 191.020), ppm = 10, mode = '-') # batch search
```

## Statistics
-----------------------
## Basic statistics

1. Find the samples names which contain the max ion intensity/peak area for each mass feature

2. Calculate coefficient of variation (CV) for each mass feature among different sample group

3. Calculate fold change (FC) for each mass feature among different sample group

4. Calculate p-value for each mass feature among different sample group

5. Calculate calculate VIP value for each mass feature among different sample group

You can use `Dostat()` function to get above statistical analysis.

``` r
# sample data
dat <- matrix(runif(2 * 300), ncol = 2, nrow = 300)
rownames(dat) <- 1:dim(dat)[1]
myGroup <- rep_len(LETTERS[1:3], 300)
# statistics
myResult <- doStat(dat, Group = myGroup)
```

## 6. PCA plot
``` r
# sample data
dat <- matrix(runif(2*300), ncol = 2, nrow = 300)
Group <- rep_len(LETTERS[1:3], 300)
# PCA
viewPCA(dat, Group = Group, interactive = T) # you can turn on/off interactive plot using interactive = T/F
```

## 7. View variations of TIC among samples

``` r
# sample data
dat <- matrix(runif(100*9), ncol = 100, nrow = 27)
myGroup <- rep_len(LETTERS[1:3], 27)
myBatch <- rep(1:3, each = 9, times = 1)
mySeq <- c(1:27)
# view TIC
viewTIC(dat, Group = myGroup, Batch = myBatch, resultBy = "Batch")
```

## 8. View volcano plot


## 9. Normalization

Normalization methods include: (1) LBME: linear baseline normalization based on mean values; (2) LBMD: linear baseline normalization based on median values; (3) PQN: probabilistic quotient normalization; (4) QT: quantile normalization; (5) TIC: total ion current normalization.

``` r
dat <- matrix(runif(100*9), ncol = 100, nrow = 10)
out <- doNormalization(dat, method = "PQN" )
```
