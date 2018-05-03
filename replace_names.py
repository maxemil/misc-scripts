#! /usr/bin/env python3
"""replace content in multiple files by patterns in another"""
import argparse
import shutil
import random

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="pattern, TAB separated")
parser.add_argument("-f", "--files", required=True, nargs='+',
                    help="files that are going to be searched and replaced")
parser.add_argument("-b", action='store_true',
                    help="switch for creating backup file with extension .bak")

args = parser.parse_args()


def replace(string, patterns):
    for k,v in patterns.items():
        string = string.replace(k, v)
    return string


def safe_replace(string, patterns):
    for p in list(patterns):
        if not any([p in key for key in patterns.keys() if not p is key]):
            string = string.replace(p, patterns[p])
            patterns.pop(p)

    while patterns:
        p = random.choice(list(patterns.keys()))
        if any([p in key for key in patterns.keys() if not p is key]):
            pass
        else:
            string = string.replace(p, patterns[p])
            patterns.pop(p)
    return string


def read_patterns(inputfile):
    patterns = {}
    for line in open(inputfile):
        line = line.strip().split('\t')
        patterns[line[0]] = line[1]
    return patterns


def need_safe_replace(patterns):
    for p in sorted(patterns.keys()):
        p_old = p
        if p_old in p:
            print("WARNING: patterns are not unique or substrings of each other! Will attempt a safe replace.")
            return True
    return False


def main(inputfile, replacementfiles, backup):
    if backup:
        for refile in replacementfiles:
            shutil.copyfile(refile, "%s.bak" % refile)
    else:
        pass

    patterns = read_patterns(inputfile)
    safe = need_safe_replace(patterns)

    for refile in replacementfiles:
        string = open(refile).read()
        if safe:
            string = safe_replace(string, patterns)
        else:
            string = replace(string, patterns)
        with open(refile, 'w') as out:
            out.write(string)


if __name__ == "__main__":
    main(args.input, args.files, args.b)
