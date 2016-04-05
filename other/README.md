# Scripts and pre-processed data

`average_shape.py`: Given a tab-delimited shape file, compute average for each position.

`correlation_analysis.html`: Jupyter notebook excerpt (Python) for performing correlating average shape profiles for different classes.

`intersect_and_split.sh`: Given a set of BED files (command line), perform intersections to generate a list of intervals that are 'unique' at specified windows. In outline, this script defines windows around each site and retrieves sites that do not overlap any other sites.

`Unique_Sites_300bp`: Abf1, Reb1, and Rap1 'unique' sites that do not overlap ChEC-defined sites in +/-300 bp windows centered at peak maxima. Files are in BED format.

