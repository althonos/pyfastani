import argparse
import glob
import json
import os
import time
import statistics
import sys
import warnings
from itertools import islice

sys.path.insert(0, os.path.realpath(os.path.join(__file__, "..", "..", "..")))

import Bio.SeqIO
import numpy
import pandas
import rich.progress
from pyfastani import Sketch

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runs", default=3, type=int)
parser.add_argument("-d", "--data", required=True)
parser.add_argument("-o", "--output", required=True)
parser.add_argument("-j", "--jobs", default=os.cpu_count() or 1, type=int)
args = parser.parse_args()

genomes = [
    list(Bio.SeqIO.parse(filename, "fasta"))
    for filename in glob.glob(os.path.join(args.data, "*.fna"))
]

with rich.progress.Progress(transient=True) as progress:
    warnings.showwarning = lambda msg, c, f, l, file=None, line=None: progress.print(msg)

    sketch = Sketch()
    task = progress.add_task(total=len(genomes), description="Sketching...")
    for genome in progress.track(genomes, task_id=task):
        sketch.add_draft(genome[0].id, [ bytes(contig.seq) for contig in genome ])
    progress.remove_task(task_id=task)

    mapper = sketch.index()

    results = dict(results=[])
    task = progress.add_task(total=len(genomes), description="Querying...")
    for genome in progress.track(genomes, task_id=task):
        contigs = [ bytes(contig.seq) for contig in genome ]
        task2 = progress.add_task(total=args.jobs, description="Threads...")
        for thread_count in progress.track(range(1, args.jobs+1), task_id=task2):
            times = []
            task3 = progress.add_task(total=args.runs, description="Repeat...")
            for run in progress.track(range(args.runs), task_id=task3):
                t1 = time.time()
                hits = mapper.query_draft(contigs, threads=thread_count)
                t2 = time.time()
                times.append(t2 - t1)
            progress.remove_task(task3)
            results["results"].append({
                "genome": genome[0].id,
                "threads": thread_count,
                "nucleotides": sum(map(len, contigs)),
                "times": times,
                "mean": statistics.mean(times),
                "stddev": statistics.stdev(times),
                "median": statistics.median(times),
                "min": min(times),
                "max": max(times),
            })
        progress.remove_task(task_id=task2)
    progress.remove_task(task_id=task)

with open(args.output, "w") as f:
    json.dump(results, f, sort_keys=True, indent=4)
