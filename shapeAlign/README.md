#shapeAlign
shapeAlign is a C++ implementation of an approach for sequence-agnostic alignment of DNA shape matrices.

[`shape_alignment_method.pdf`](shape_alignment_method.pdf) contains an outline of how the shape alignment approach works.

`test/` contains a Perl script for generating simulated data for alignment.

**Compilation & Dependencies**

Typing `make` should compile shapeAlign. Note that [GSL](http://www.gnu.org/software/gsl/) is required (GSL matrices and the GSL BLAS interface are used extensively).

**Usage**

`./shapeAlign` without any arguments should display a help message.

General usage:

```
Usage:   shapeAlign [OPTIONS] -n <site name file> -f <shape file 1> ... <shape file n> 

 *Note: shape file list must be the last entered command line argument.

 -n      Single-column file containing names/accessions of sites.
 -f      Space-delimited list of shape files (multi-column, tab-
         delimited file). By default, the first three columns and
         the last two columns are ignored. The ordering of sites
         in these files is assumed to be the same as in the site
         names file.

 Options: 
 -min    Minimum shift. (Default: -25 bp).
 -max    Maximum shift. (Default: 25 bp).
 -start  Start of alignment window defined with respect to the
         start of the shape window, which is defined as zero.
         (Default: 0).
 -len    Length of the alignment window. (Default: length of 
         shape window).
```