#! /bin/bash


# User: fcruz
# Date: 2022-06-09

# Brief description/notes:
# Pass input via standard input (gap file should preserve original contig names, use fasta-stats-py)

# Dependencies


gawk '{print $1"\t"$2"\t"$3"\t100"}'
