import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

with open(args.input, "r") as src:
    source = src.read()

# make everything public!
source = re.sub("private:", "public:", source)
# patch some classes with private members at the
# beginning of their declaration (skch::Sketch)
source = re.sub(r"class (.+\s*)\{(\s*)", r"class \1 {\npublic:\n\2", source)
# patch the contructors of skch::Map and skch::Seq so that they don't
# do anything
source = re.sub(r"(Sketch|Map)\(([^)]*)\)([^}]+)\{[^}]*}", r"\1(\2)\3 {}", source)

os.makedirs(os.path.dirname(args.output), exist_ok=True)
with open(args.output, "w") as dst:
    dst.write(source)