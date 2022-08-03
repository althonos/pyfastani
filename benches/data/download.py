#!/bin/sh

import os
import gzip
import shutil
import urllib.request

import rich.progress

SAMPLES = [
    "76859.SAMN03263149",
    "76856.SAMN03263147",
    "562982.SAMN02595332",
    "154046.SAMEA3545258",
    "322095.SAMN03854412",
    "596327.SAMN00002220",
    "28131.SAMD00034934",
    "1123263.SAMN02441153",
    "33039.SAMEA3545241",
    "1512.SAMEA3545330",
    "926561.SAMN02261313",
    "1121279.SAMN02745887",
    "1216006.SAMEA4552037",
    "715226.SAMN02469919",
    "1121935.SAMN02440417",
    "1210082.SAMD00040656",
    "1189612.SAMN02469969",
    "518637.SAMN00008825",
    "585506.SAMN00139194",
    "82977.SAMN02928654",
    "1267533.SAMN02440487",
    "373.SAMN04223557",
    "298701.SAMN02471586",
    "1852381.SAMEA4555866",
    "1341181.SAMN02470908",
    "1194088.SAMN02470121",
    "198092.SAMN02745194",
    "1561204.SAMN03107784",
    "570952.SAMN02441250",
    "419481.SAMN05216233",
    "888845.SAMN05017598",
    "1254432.SAMN02603858",
    "1909395.SAMN05912833",
    "448385.SAMEA3138271",
    "1343740.SAMN02604192",
    "48.SAMN03704078",
    "1391654.SAMN03951129",
    "1912.SAMN05935554",
    "576784.SAMEA2519540",
    "749414.SAMN02603683",
    "722472.SAMN05444171",
    "1173028.SAMN02261253",
    "1134687.SAMN02463898",
    "398579.SAMN02598386",
    "54388.SAMEA2645827",
    "1773.SAMN06010452",
    "1281227.SAMN01885939",
    "645465.SAMN02595349",
    "36809.SAMN04572937",
    "649639.SAMN00007429",
]

data_folder = os.path.dirname(os.path.realpath(__file__))
for sample in rich.progress.track(SAMPLES, description="Downloading..."):
    tax_id = sample.split(".")[0]
    url = "https://progenomes.embl.de/dumpSequence.cgi?p={}&t=c&a={}".format(sample, tax_id)
    with urllib.request.urlopen(url) as res:
        with gzip.open(res) as src:
            with open(os.path.join(data_folder, "{}.fna".format(sample)), "wb") as dst:
                shutil.copyfileobj(src, dst)
