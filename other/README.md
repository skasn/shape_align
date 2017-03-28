# Scripts and pre-processed data

Instructions for use:

1. Get BED files (available in [`sites/bed`](https://github.com/sivakasinathan/shape_align/tree/master/other/sites/bed)) and FASTA files (available in [`sites/fasta`](https://github.com/sivakasinathan/shape_align/tree/master/other/sites/fasta)). The number of sites for each window size used the analysis is available in [`sites/num_sites.txt`](https://github.com/sivakasinathan/shape_align/blob/master/other/sites/num_sites.txt).

1. Get sequences from the FASTA files (referenced above) corresponding to the sites in the BED files (also referenced above). Create a separate FASTA file for each of the desired analysis classes (i.e., “Fast” and “Slow” files for each TF). A desired number of random sites can be sampled from the files provided.

1. Get shape features for sites using [DNAshapeR](https://www.ncbi.nlm.nih.gov/pubmed/26668005) from the Rohs lab. [Detailed instructions](http://rohslab.cmb.usc.edu/Documents/DNAshapeR_document.pdf) for using DNAshapeR with typical use-case vignettes are available.

1. Convert DNAshapeR comma-separated output to tab-delimited format. The first column should be the site name or other unique identifier and subsequent columns should contain the values for the shape parameters.

1. Install dependencies for shapeAlign following the instructions from the distributors of the following packages: [OpenMP](http://openmp.org/) and [GSL](http://www.gnu.org/software/gsl/).

1. Compile shapeAlign: Once the dependencies are installed, download the [shapeAlign source code](https://github.com/sivakasinathan/shape_align/tree/master/shapeAlign), navigate to the directory containing the makefile packaged with shapeAlign and type ‘make.’ Issues during compilation with dependencies can be resolved by modifying the makefile by including the appropriate paths to the GSL/OpenMP dependencies or the version of the g++ compiler indicated in the makefile.

1. Generate ‘names’ file for shapeAlign – this is simply the first column (containing the site names/unique identifiers) of one of the shape parameter files. It does not matter which file is used to generate this file since the shape parameters in the files must be ordered in the same way.

1.  Run shapeAlign using the following parameters:

    `./shapeAlign -start 50 -end 141 -min -25 -max 25 -istart 97 -iend 102 –n [path to names file] –f [shape file 1] … [shape file 4]`

1. Determine averages for shape features ignoring ‘NA’ values and correlate average profiles using, e.g., Excel, R, or Python/Pandas.

1. Correlations of average profiles can be performed using the code provided in `correlation_analysis/correlation_analysis.html`. Heatmaps with the correlation values included are in `correlation_analysis/*.png`.

---

Some useful scripts:

`average_shape.py`: Given a tab-delimited shape file, compute average for each position.

`intersect_and_split.sh`: Given a set of BED files (command line), perform intersections to generate a list of intervals that are 'unique' at specified windows. In outline, this script defines windows around each site and retrieves sites that do not overlap any other sites.

