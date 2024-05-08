# Data Processing Files

These files are mostly meant to analyze and process data as pre-processing step to the model. They are run once on data, and then data is used for the model. 

## Description of Relevant Files:

### ReadStatistics 
- calculates read-level stats like average read length; note: requires files in FASTQ format
- python3 ReadStatistics_shree.py <path to FASTQ file>

### AlignmentStatistics 
- calculates alignment-level stats like number/percentage of ambiguous reads; note: requires files in BAM/SAM format
- python3 AlignmentStatistics_shree.py <path to BAM/SAM file>

### ParseTxt 
- takes alignments in TXT format and creates a JSON file with the same format as that created from a BAM file; requires file in TXT format; path to JSON file (including file name) should be specified 
- python3 ParseTxt.py <path to text file> <path where JSON file shoule be created including JSON file name>
- example: python3 ParseTxt.py PacBio/PacBio/day0_rep2/collapsed_transcripts.read_stat.txt PacBio/PacBio/day0_rep2/day0_rep2.json

### Downsample
- this is meant to be run only once
- file paths are already specified in file
- downsamples alignment files according to algorithm specified in report
- python3 downsample.py
