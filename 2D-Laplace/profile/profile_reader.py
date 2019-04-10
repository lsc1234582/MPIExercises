import argparse
import numpy as np
import os
import sys
import pandas as pd
import re
import glob
import shutil as sh
import json

def curdir():
    return os.path.abspath(os.curdir)

def discover_profile_Cases(path):
    os.chdir(path)
    dir_stack = []
    case_freq = {}
    cases = []
    for i, p in enumerate(glob.iglob("**/profilecase_*", recursive=True)):
        case_dir = os.path.abspath(p)
        dir_stack.append(curdir())
        os.chdir(case_dir)

        # Parse a case
        try:
            with open("profile.meta/meta.json", "r") as fh:
                case_meta = json.load(fh)
        except Exception as e:
            os.chdir(dir_stack.pop())
            continue

        case_entry = case_meta.copy()
        del case_entry["CommandTimeout"]
        del case_entry["Comments"]
        with open("profile.json", "r") as fh:
            case_profile = json.load(fh)
            case_entry["Machine"] = case_profile["info"]["machine"]
            del case_profile

        #print(case_dir)
        case_entry["CasePath"] = case_dir
        cases.append(case_entry)
        #print(case_entry)
        case_label = tuple([(k, v) for k, v in case_entry.items()])
        freq = case_freq.setdefault(case_label, 0)
        case_freq[case_label] += 1
        os.chdir(dir_stack.pop())

    cases = pd.DataFrame(cases)
    return cases

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Laplace solver profile reader")
    parser.add_argument("-r", "--root_path", dest="root_path", metavar="Profile root path", action="store", type=str,
            help=("Path to a root directory which contains all the profile cases"), required=True)
    parser.add_argument("-o", "--output_name", dest="output_name", metavar="Output csv name", action="store", type=str,
            help=("Name of the output csv file"), required=True)
    args = parser.parse_args()
    cases = discover_profile_Cases(args.root_path)
    cases.to_csv(args.output_name)
