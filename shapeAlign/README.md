#shapeAlign
shapeAlign is a C++ implementation of an approach for sequence-independent alignment of DNA shape matrices.

[`shape_alignment_method.pdf`](shape_alignment_method.pdf) contains an outline of how the shape alignment approach works.

`test/` contains a Perl script for generating simulated data for alignment.

**Compilation & Dependencies**

Typing `make` should compile shapeAlign. *Note* that [GSL](http://www.gnu.org/software/gsl/) is required (GSL matrices and the GSL BLAS interface are used extensively).
[OpenMP](http://openmp.org/wp/) is used for parallelization of centroid finding.

**Usage**

`./shapeAlign` without any arguments displays the following help message:

```
Sequence-independent alignment of DNA shape

Usage:   shapeAlign [OPTIONS] -n <site name file> -f <shape file 1> ... <shape file n>

 *Note: shape file list must be the last entered command line argument.

 -n      Single-column file containing names/accessions of sites.
 -f      Space-delimited list of shape files (multi-column, tab-
         delimited file). By default, the first three columns and
         the last two columns are ignored. The ordering of sites
         in these files is assumed to be the same as in the site
         names file.

 Options:
 -entr   Minimum number of valid bases in overlap. (Default: 10).
 -min    Minimum shift. (Default: -25 bp).
 -max    Maximum shift. (Default: 25 bp).
 -start  Start of alignment window defined with respect to the
         start of the shape window, which is defined as zero.
         (Note: If -start is used, -end must be defined).
 -end    End of alignment window defined with respect to the
         start of the shape window, which is defined as zero.
         (Note: If -end is used, -start must be defined).
 -istart  Start of positions to ignore defined with respect to
         start of the shape window, which is defined as zero.
         (Note: If -istart is used, -iend must be defined).
 -iend    End of positions to ignore defined with respect to
         start of the shape window, which is defined as zero.
         (Note: If -iend is used, -istart must be defined).
```

**Testing/Validation**

A mock dataset of simulated vectors was used to validate the alignment approach. Scripts used to validate the approach are included in `test/`.

**Analysis parameters**

For analyses of Abf1, Reb1, and Rap1 ChEC data, the following parameters were used:
`-start 50 -end 141 -min -25 -max 25 -istart 97 -iend 102`
