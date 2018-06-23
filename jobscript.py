#!python

#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# do something useful with the threads
threads = job_properties[threads]

# access property defined in the cluster configuration file (Snakemake >=3.6.0)
job_properties["cluster"]["time"]

if not threads:
    threads = 1

os.system(
    "qsub -v PATH='/homes/22/jc220896/miniconda3/envs/virusMAP/bin:$PATH' \
    -d /homes/22/jc220896/virusMAP/ \
    -o /homes/22/jc220896/virusMAP/qsublogs/ \
    -e /homes/22/jc220896/virusMAP/qsublogs/ \
    -l pmem=2GB \
    -l nodes=1:ppn={threads} \
    {script}".format(
        threads=threads, script=jobscript)
    )

snakemake -p \
--latency-wait 60 \
--cluster  \
-j 10
