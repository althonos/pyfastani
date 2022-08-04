import argparse
import itertools
import json
import os
import random
import re
import math

import numpy
import matplotlib.pyplot as plt
import scipy.stats
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial
from palettable.cartocolors.qualitative import Prism_10

random.seed(0)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output")
parser.add_argument("-s", "--show", action="store_true")
args = parser.parse_args()

def exp_decay(x, a, b, c):
    return a*numpy.exp(-b*x) + c

with open(args.input) as f:
    data = json.load(f)

# fig = plt.figure(1)
# ax = fig.add_subplot(111, projection='3d')
#
# X = [ bench["nucleotides"] / 1e6 for bench in data["results"] ]
# Y = [ bench["threads"] for bench in data["results"] ]
# Z = [ bench["mean"] for bench in data["results"] ]
# c = [ Bold_9.hex_colors[ X.index(bench["nucleotides"] / 1e6) % 9 ] for i, bench in enumerate(data["results"]) ]
#
# ax.scatter(X, Y, Z, c=c)
# ax.set_xlabel("Genome size (Mbp)")
# ax.set_ylabel("Number of Threads")
# ax.set_zlabel("Time (s)")

plt.figure(1, figsize=(12, 6))

plt.subplot(1, 2, 1)
data["results"].sort(key=lambda r: (r["threads"], r["nucleotides"]))
for color, (threads, group) in zip(
    itertools.cycle(Prism_10.hex_colors),
    itertools.groupby(data["results"], key=lambda r: r["threads"])
):
    group = list(group)
    X = numpy.array([r["nucleotides"] / 1e6 for r in group])
    Y = numpy.array([r["mean"] for r in group])

    # p = Polynomial.fit(X, Y, 2)
    # pX = numpy.linspace(0, max(r["sequences"] for r in group), 1000)
    # reg = scipy.stats.linregress(X, Y)
    # plt.plot([ 0, max(X) ], [ reg.intercept, reg.slope*max(X) + reg.intercept ], color=color, linestyle="--", marker="")
    # ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group]
    plt.scatter(X, Y, marker="+", color=color, label=f"Threads={threads}")
    # plt.plot(pX, p(pX), color=color, linestyle="--")

plt.legend()
plt.xlabel("Genome size (Mbp)")
plt.ylabel("Query Time (s)")

plt.subplot(1, 2, 2)

values = list({ r["nucleotides"] for r in data["results"] })
choices = { random.choice(values) for _ in range(5) }
data["results"] = [ r for r in data["results"] if r["nucleotides"] in choices ]

data["results"].sort(key=lambda r: (r["nucleotides"], r["threads"]))
for color, (nucleotides, group) in zip(
    itertools.cycle(Prism_10.hex_colors),
    itertools.groupby(data["results"], key=lambda r: r["nucleotides"])
):
    group = list(group)
    X = numpy.array([r["threads"] for r in group if r["mean"]])
    Y = numpy.array([r["mean"] for r in group])

    popt, pcov = curve_fit(exp_decay, X, Y)
    pX = numpy.linspace(1, max(r["threads"] for r in group), 100)

    name = ".".join(group[0]['genome'].split(".")[:2])
    plt.scatter(X, Y, marker="+", color=color, label=f"{name} ({nucleotides/1e6:.1f} Mbp)")
    plt.plot(pX, exp_decay(pX, *popt), color=color, linestyle="--", linewidth=0.3)

plt.legend()
plt.xlabel("Number of Threads")
plt.ylabel("Query Time (s)")

plt.tight_layout()
output = args.output or args.input.replace(".json", ".svg")
plt.savefig(output, transparent=True)
if args.show:
    plt.show()
