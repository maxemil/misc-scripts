#! /usr/bin/env python3
"""replace content in multiple files by patterns in another"""
import argparse
import shutil
import random
import copy

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="pattern, TAB separated")
parser.add_argument("-f", "--files", required=True, nargs='+',
                    help="files that are going to be searched and replaced")
parser.add_argument("-b", action='store_true',
                    help="switch for creating backup file with extension .bak")
parser.add_argument("-r", "--reverse", action='store_true',
                    help="patterns are reversed (new\told) instead of (old\tnew)")

args = parser.parse_args()


def replace(filecontent, patterns):
    for k,v in patterns.items():
        filecontent = filecontent.replace(k, v)
    return filecontent


def safe_replace(filecontent, patterns):
    for p in list(patterns):
        if not any([p in key for key in patterns.keys() if not p is key]):
            filecontent = filecontent.replace(p, patterns[p])
            patterns.pop(p)

    while patterns:
        p = random.choice(list(patterns.keys()))
        if any([p in key for key in patterns.keys() if not p is key]):
            pass
        else:
            filecontent = filecontent.replace(p, patterns[p])
            patterns.pop(p)
    return filecontent


def read_patterns(inputfile, reverse):
    patterns = {}
    for line in open(inputfile):
        line = line.strip().split('\t')
        if reverse:
            patterns[line[1]] = line[0]
        else:
            patterns[line[0]] = line[1]
    return patterns


def need_safe_replace(patterns):
    p_old = ""
    for p in sorted(patterns.keys()):
        if not p_old:
            p_old = p
            continue
        if p_old in p:
            print("WARNING: patterns are not unique or substrings of each other! Will attempt a safe replace.")
            print(p_old, p)
            return True
        p_old = p
    return False


def main(inputfile, replacementfiles, backup, reverse):
    if backup:
        for refile in replacementfiles:
            shutil.copyfile(refile, "%s.bak" % refile)
    else:
        pass

    patterns = read_patterns(inputfile, reverse)
    safe = need_safe_replace(patterns)

    for refile in replacementfiles:
        filecontent = open(refile).read()
        if safe:
            filecontent = safe_replace(filecontent, copy.deepcopy(patterns))
        else:
            filecontent = replace(filecontent, patterns)
        with open(refile, 'w') as out:
            out.write(filecontent)


if __name__ == "__main__":
    main(args.input, args.files, args.b, args.reverse)
