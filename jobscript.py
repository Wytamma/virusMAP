#!python

#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)


# do something useful with the threads
if threads in job_properties:
    threads = job_properties[threads]
else:
    threads = 1

# access property defined in the cluster configuration file (Snakemake >=3.6.0)
#job_properties["cluster"]["time"]


os.system(
    """qsub -v PATH='/homes/22/jc220896/miniconda3/envs/virusMAP/bin:$PATH' \
    -l pmem=2GB \
    -l nodes=1:ppn={{threads}} \
    {{jobscript}}""".format(jobscript=jobscript, threads=threads)
