import argparse
import numpy as np
import os
import sys
import pandas as pd
import re
import glob
import shutil as sh
import matplotlib.pyplot as plt
import json
import seaborn as sns
import tarfile


# Solver category order mapping: Axis type -> key func
solver_map = {k:v + 1 for v, k in enumerate(["Serial", "ParRow", "Par"])}
func_map = {k:v + 1 for v, k in enumerate(["0", "1", "2", "3"])}

axis_map_map = {
    "Solver": lambda v: solver_map[v],
    "Func": lambda v: func_map[v],
    "NRow": lambda v: v,
    "NCol": lambda v: v,
    "NPX" : lambda v: v,
    "NPY" : lambda v: v,
}

# Metrics to ylabel mapping
metrics_2_ylabel = {
    "RunTime": "RunTime / ms",
    "Compute+MPI+IO": "Compute+MPI+IO percentage / %"
}

def curdir():
    return os.path.abspath(os.curdir)

def unique_no_sort(a):
    indexes = np.unique(a, return_index=True)[1]
    return [a[index] for index in sorted(indexes)]

def get_case_parameters(cases, legacy=False):
    parameters_exclude=["EndDateTime","StartDateTime"]
    if legacy:
        parameters_exclude.append("CasePath")

    return {param: list(np.unique(cases.loc[:, param])) for param in cases.columns.drop(parameters_exclude)}

def transform_with_labels(df_tuples, axes, metrics):
    """Genereate and add axis labels, x-axis tick labels, category keys/index to the list of tuples,
    and sort the tuple according to the category key.
    """
    def entry_key(val):
        return val[0] * val[0][1]

    if len(axes) == 2:
        print("HAAAAAAAAAA")
        df_tuples_with_label = [[i, axis_map_map[axes[0]](i), v, j] for (i, j), v in df_tuples]
        xlabel = axes[0]

    elif len(axes) == 1:
        if "^" in axes[0]:
            # Composite axes
            sub_axes = axes[0].split("^")
            if len(sub_axes) != 2:
                raise Exception("Number of sub-axes in a composite axis has to be 2")
            df_tuples_with_label = [[str(i) + "^" + str(j), i * j, v] for (i, j), v in df_tuples]
            # Normlize category indices
            norm_is = normalize_ns([i for k, i, v in df_tuples_with_label])
            df_tuples_with_label = [[k, norm_is[i], v] for i, (k, _, v) in enumerate(df_tuples_with_label)]
            # Combine xtick labels if landing on same index
            index_label_dict = {}
            for ind, (k, i, v) in enumerate(df_tuples_with_label):
                index_label_dict.setdefault(i, [])
                index_label_dict[i].append((ind, k))
            for k, v in index_label_dict.items():
                indices, keys = zip(*v)
                combined_key = ",".join(np.unique(keys))
                for ind in indices:
                    df_tuples_with_label[ind][0] = combined_key
            xlabel = str(sub_axes[0]) + "x" + str(sub_axes[1])
        else:
            # Singular axis
            df_tuples_with_label = [[k, axis_map_map[axes[0]](k), v] for k, v in df_tuples]
            xlabel = axes[0]
    else:
        raise Exception("Does not support more than 2 axes")


    # Sort categories according to category index
    df_tuples_with_label.sort(key=(lambda v: v[1]))

    return df_tuples_with_label, xlabel, metrics_2_ylabel[metrics]

def is_integral(t):
    return t in [int, np.int64, np.int32]

def select_cases(cases, sel_dict, composite_axes, metrics, legacy_load=False, copy_cases=False):
    #group_by_axes = axes[0].split("^") if len(axes) == 1 and "^" in axes[0] else axes
    if copy_cases:
        plot_df = cases.copy()
    else:
        plot_df = cases

    if len(sel_dict) != 0:
        sel = list(map(list, zip(*list(sel_dict.items()))))
        sel = np.all(cases.loc[:, sel[0]] == sel[1], axis=1)

        plot_df = plot_df.loc[sel]

    # Composite axes : composite axis name -> list(singular axes names)
    c_axes = {}
    for c_ax in composite_axes:
        assert("^" in c_ax)
        c_ax = c_ax.split(":")
        if len(c_ax) == 2:
            c_ax_name, c_ax = c_ax
        elif len(c_ax) == 1:
            c_ax_name = c_ax[0]
            c_ax = c_ax[0]
        else:
            assert(False)
        axes = c_ax.split("^")
        assert(len(axes) == 2)
        c_axes[c_ax_name] = axes

    for i, row in plot_df.iterrows():
        # Add composite axis
        for c_ax_name, axes in c_axes.items():
            v1 = plot_df.loc[i, axes[0]]
            v2 = plot_df.loc[i, axes[1]]
            if is_integral(type(v1)) and is_integral(type(v2)):
                plot_df.loc[i, c_ax_name] = v1 * v2
            else:
                #print("hahaha")
                #print(row)
                #print(c_ax_name)
                #print("hahaha")
                plot_df.loc[i, c_ax_name] = str(v1) + "^" + str(v2)
                #print(plot_df.loc[i, c_ax_name])
                #print("hahaha2")

        if legacy_load:
            try:
                with open(os.path.join(row["CasePath"], "profile.json"), "r") as fh:
                    profile_json = json.load(fh)
                add_metrics(metrics, plot_df, i, profile_json)
                del profile_json
            except Exception as e:
                timeout_init = os.path.exists(os.path.join(row["CasePath"], "TIMEOUT_INIT"))
                timeout_solve = os.path.exists(os.path.join(row["CasePath"], "TIMEOUT_SOLVE"))
                if not timeout_init and not timeout_solve:
                    print("ERROR: {} misses profile.json".format(row["CasePath"]))
                    raise e
        # Make composite categorical
        for c_ax_name, axes in c_axes.items():
            c_ax_col = plot_df.loc[:, c_ax_name]
            if c_ax_col.dtype != np.int and c_ax_col.dtype != np.float:
                plot_df.loc[:, c_ax_name] = c_ax_col.astype("category")

    if "Speedup" in metrics:
        max_runtime = np.max(plot_df.loc[:, "RunTime"])
        plot_df = plot_df.assign(Speedup=lambda x: max_runtime / x["RunTime"])
    if "RunTimeNorm" in metrics:
        min_runtime = np.min(plot_df["RunTime"])
        plot_df = plot_df.assign(RunTimeNorm=lambda x: x["RunTime"] / min_runtime)
        plot_df["RunTimeNorm"] = plot_df["RunTime"] / min_runtime

    #group_by_axes.append(metrics)
    #print(group_by_axes)
    return plot_df

def generate_plot_dict(cases, sel_dict, axes, metrics):
    group_by_axes = axes[0].split("^") if len(axes) == 1 and "^" in axes[0] else axes

    sel = list(map(list, zip(*list(sel_dict.items()))))
    #print(sel)
    sel = np.all(cases.loc[:, sel[0]] == sel[1], axis=1)

    plot_df = cases.loc[sel].copy()
    #print(plot_df)
    for i, row in plot_df.iterrows():
        #print(row["CasePath"])
        with open(os.path.join(row["CasePath"], "profile.json"), "r") as fh:
            profile_json = json.load(fh)
        # Add value
        if metrics == "RunTime":
            plot_df.loc[i, metrics] = profile_json["info"]["runtime"]
        #plot_df
    #plot_df["a"] = 2
    #cases.loc[pd.notna()["Func"], "Func"]
    plot_df.describe()
    plot_df_dict = {k: list(v.loc[:, metrics]) for k, v in plot_df.groupby(group_by_axes)}
    plot_df_tuple= [(k, v) for k in plot_df_dict.keys() for v in plot_df_dict[k]]

    plot_df_tuple_with_label, xlabel, ylabel = transform_with_labels(plot_df_tuple, axes, metrics)
    plot_dict = {
        "xlabel": xlabel,
        "ylabel": ylabel,
        "df": plot_df_tuple_with_label,
        "sel": sel_dict
    }

    return plot_dict

def plot(plot_dict):
    plot_df_series = list(zip(*(plot_dict["df"])))
    print((unique_no_sort(plot_df_series[1]), unique_no_sort(plot_df_series[0])))
    plt.xticks(unique_no_sort(plot_df_series[1]), unique_no_sort(plot_df_series[0]), rotation='vertical')
    #plt.xticks(plot_df_series[1], plot_df_series[0])

    plt.xlabel(plot_dict["xlabel"])
    plt.ylabel(plot_dict["ylabel"])
    #print((plot_df_series[1], plot_df_series[2]))
    if len(plot_df_series) == 3:
        # Single-axis plot
        sns.scatterplot(plot_df_series[1], plot_df_series[2], alpha=0.9)
    elif len(plot_df_series) == 4:
        # Double-axis plot
        sns.scatterplot(plot_df_series[1], plot_df_series[2], hue=plot_df_series[3], alpha=0.9)
    #print(plot_dict["sel"])
#plot_df.loc[:, axes + [metrics]].groupby(axes).agg([np.mean, np.std])
#plot_df.loc[:,["NPX", "NPY", "RunTime"]]

def abspath(p):
    if os.path.isabs(p):
        return p
    return os.path.join(curdir(), p)

def generate_profile_cases(root, uncompress_folder=".profiletmp", delete_uncompress_folder=False, uncompress_member_func=None):
    """TODO:
    """
    dir_stack = []
    cases = []
    dir_stack.append(curdir())
    os.chdir(root)
    if os.path.exists(uncompress_folder):
        raise Exception("Uncompress folder {} already exists".format(uncompress_folder))
    os.makedirs(uncompress_folder)

    def is_tar_gz_file(f):
        return os.path.isfile(f) and tarfile.is_tarfile(f) and "gz" in os.path.splitext(f)[1]

    for p in glob.iglob("**/profile_*", recursive=True):
        print("INFO: profile folder/file {} detected".format(p))
        if p.split(os.path.sep)[0] == uncompress_folder:
            # Make sure do not re-glob uncompressed folders
            continue
        elif os.path.isdir(p):
            profile_dir = abspath(p)
        elif is_tar_gz_file(p):
            p_without_ext = os.path.splitext(p)[0]
            p_without_ext = os.path.splitext(p_without_ext)[0]
            tmp_p_dir = os.path.join(uncompress_folder, p_without_ext)
            os.makedirs(tmp_p_dir)
            with tarfile.open(p) as t:
                if uncompress_member_func is None:
                    t.extractall(path=uncompress_folder)
                else:
                    t.extractall(members=uncompress_member_func(t), path=uncompress_folder)
            profile_dir = abspath(tmp_p_dir)
        else:
            continue

        dir_stack.append(curdir())
        os.chdir(profile_dir)
        yield from glob.iglob("**/profilecase_*", recursive=True)

        if delete_uncompress_folder and is_tar_gz_file(p):
            p_without_ext = os.path.splitext(p)[0]
            p_without_ext = os.path.splitext(p_without_ext)[0]
            tmp_p_dir = os.path.join(uncompress_folder, p)
            sh.rmtree(tmp_p_dir, ignore_errors=True)
        os.chdir(dir_stack.pop())
    if delete_uncompress_folder:
        sh.rmtree(uncompress_folder, ignore_errors=True)
    return pd.DataFrame(cases)


def load_cases(root, metrics):
    """Discover all profile cases under a root dir. The profile cases can be in a folder or a tar.gz compressed archive
    Return a DataFrame object containing all the params combination and a path to each case.
    """
    dir_stack = []
    cases = []

    def profile_and_meta_json_file(members):
        for tarinfo in members:
            if "profile_" in tarinfo.name and\
                    "profilecase_" in tarinfo.name and\
                    ("profile.meta/meta.json" in tarinfo.name or\
                     "profile.json" in tarinfo.name):
                yield tarinfo

    for case_dir in generate_profile_cases(root, ".profiletmp", True, profile_and_meta_json_file):
        dir_stack.append(curdir())
        os.chdir(case_dir)
        # Parse a case
        try:
            with open("profile.meta/meta.json", "r") as fh:
                case_entry = json.load(fh)
        except Exception as e:
            print("ERROR: {} misses meta.json".format(case_dir))
            os.chdir(dir_stack.pop())
            raise e
        del case_entry["SolverTimeout"]
        del case_entry["Comments"]

        timeout_init = False
        timeout_solve = False
        try:
            with open("profile.json", "r") as fh:
                profile_json = json.load(fh)
        except Exception as e:
            timeout_init = os.path.exists("TIMEOUT_INIT")
            timeout_solve = os.path.exists("TIMEOUT_SOLVE")
            if not timeout_init and not timeout_solve:
                print("WARNING: {} misses profile.json".format(case_dir))
                os.chdir(dir_stack.pop())
                continue
        case_entry["Machine"] = profile_json["info"]["machine"]
        case_entry["TimeoutInit"] = timeout_init
        case_entry["TimeoutSolve"] = timeout_solve
        # Add metrics
        add_metrics_dict(metrics, case_entry, profile_json)
        #print(case_dir)
        cases.append(case_entry)
        del profile_json
        #print(case_entry)
        os.chdir(dir_stack.pop())
        # TODO: Remove current case after being done with it (if uncompressed)?
        #sh.rmtree(case_dir, ignore_errors=True)
    return pd.DataFrame(cases)

def add_metrics(metrics, cases, ind, profile_json):
    for m in metrics:
        if m == "RunTime":
            cases.loc[ind, m] = profile_json["info"]["runtime"]
        elif m == "ComputeTime":
            cases.loc[ind, m] = np.mean(profile_json["samples"]["activity"]["main_thread"]["normal_compute"])
        elif m == "MPITime":
            cases.loc[ind, m] =\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["point_to_point_mpi"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["point_to_point_mpi_non_main_thread"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["point_to_point_mpi_openmp"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["collective_mpi"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["collective_mpi_non_main_thread"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["collective_mpi_openmp"])
        elif m == "IOTime":
            cases.loc[ind, m] =\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["io_reads"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["io_reads_openmp"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["io_writes"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["io_writes_openmp"])
        elif m == "SampleSize":
            cases.loc[ind, m] = len(profile_json["samples"]["activity"]["main_thread"]["io_reads"])

#help(glob)
def add_metrics_dict(metrics, case_entry, profile_json):
    for m in metrics:
        if m == "RunTime":
            case_entry[m] = profile_json["info"]["runtime"]
        elif m == "ComputeTime":
            case_entry[m] = np.mean(profile_json["samples"]["activity"]["main_thread"]["normal_compute"])
        elif m == "MPITime":
            case_entry[m] =\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["point_to_point_mpi"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["point_to_point_mpi_non_main_thread"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["point_to_point_mpi_openmp"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["collective_mpi"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["collective_mpi_non_main_thread"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["collective_mpi_openmp"])
        elif m == "IOTime":
            case_entry[m] =\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["io_reads"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["io_reads_openmp"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["io_writes"]) +\
                    np.mean(profile_json["samples"]["activity"]["main_thread"]["io_writes_openmp"])
        elif m == "SampleSize":
            case_entry[m] = len(profile_json["samples"]["activity"]["main_thread"]["io_reads"])

#help(glob)

def preprocess_cases(cases, legacy=False):
    # Pre-processing
    # Fill in missing values (NPX, NPY)
    if legacy:
        cases.loc[(cases.loc[:, "NPX"]) == "", "NPX"] = 1
        cases.loc[(cases.loc[:, "NPY"]) == "", "NPY"] = 1
        cases.loc[:, "machine"] = cases.loc[:, "machine"].fillna("e121008-lin")
        cases = cases.assign(SolverTolerance=0.00000001,
                XMin=1.0,
                XMax=-1.0,
                YMin=1.0,
                YMax=-1.0)
    # Transform type
    if legacy:
        cases.loc[:, ["NCol", "NRow", "NPX", "NPY"]] = cases.loc[:, ["NCol", "NRow", "NPX", "NPY"]].applymap(int)
    # Make categories
    categorical_axes = ["AdditionalTags", "BuildCommand", "Func", "Machine", "Solver", "SolverLaunchCommand", "SourceVersionTag",
            "TimeoutInit", "TimeoutSolve"]
    cases.loc[:, categorical_axes] = cases.loc[:, categorical_axes].astype("category")
    return cases


def gen_fullfac_from_basefac(bfs):
    """ Generate a full list of unique factors from the factorisation dict that contains the frequencies of base(prime)factors
    """
    bfs.sort()
    return list(do_gen_fullfac_from_basefac(bfs, 0, True, 1, 0, np.ceil(len(bfs) / 2.0) - 1, np.prod(bfs)))


def do_gen_fullfac_from_basefac(bfs, curr_loc, parent_leftmost, fac, depth, max_depth, original_product):
    """ A somewhat n-ary-tree approach, recursive.
    bfs = [p1,       p1, p1, ... p2, ..., pn, ..., pn]
    """
    next_locs = np.unique(bfs, return_index=True)[1]
    for i, next_loc in enumerate(next_locs):
        # Next becomes current
        new_fac = fac * bfs[next_loc]
        yield (new_fac, int(original_product / new_fac))
        # current becomes parrent
        new_bfs = bfs[next_loc:].copy()
        new_bfs.pop(0)
        parent_leftmost &= (i == 0)
        if (depth < max_depth - 1 or (depth == max_depth - 1 and parent_leftmost)):
            yield from do_gen_fullfac_from_basefac(new_bfs, next_loc, parent_leftmost,
                                                                   new_fac, depth + 1, max_depth, original_product)

def do_gen_fullfac_from_basefac_iter(bfs):
    """ A somewhat n-ary-tree approach, iterative.
    bfs = [p1, p1, p1, ... p2, ..., pn, ..., pn]
    """
    pass

def fac_brute(n):
    """ Factorise a number n into a list of factors/divisors
    Brute-force solution
    """
    return [f for f in range(1, n + 1) if n % f == 0]

def fac_sym(n):
    """ Factorise a number n into a list of factors/divisors
    Brute-force solution; takes advantages of symmetry
    """
    return [(f1, int(n / f1)) for f1 in range(1, int(np.sqrt(n)) + 1) if n % f1 == 0]

def fac_prime(n, ps):
    """ Factorise a number n into a list of factors/divisors
    Use a list of prime numbers to facilitate factorisation
    """
    cur_n = n
    fs = [(1, n)]
    while cur_n != 1:
        print(cur_n)
        for p in ps:
            f = cur_n % p
            cum_p = p
            while f == 0:
                fs.append((cum_p, int(n / cum_p)))
                cur_n //= p
                cum_p *= p
                f = cur_n % p
    return fs

def gcd_brute(ns):
    """Find greatest common divisor of a list of (unsorted) numbers
    Brute-force solution
    [1, 2, 3, 90, ]
    [3, 6, 19,]
    [2, 9, 13,]
    """
    # Generate a list of list of factors for each number
    fss = []
    for n in ns:
        fs = list(zip(*fac_sym(n)))
        fs_list = list(fs[0])
        fs_list.extend(fs[1])
        fs_list.sort(reverse=True)
        fss.append(fs_list)

    # Find the gcd from fss
    empty_fs_list = 0
    while empty_fs_list < len(fss):
        max_cd = 0
        max_cd_inds = []
        for i, fs in enumerate(fss):
            if len(fs) == 0:
                continue
            f = fs[0]
            if max_cd < f:
                max_cd = f
                max_cd_inds.clear()
                max_cd_inds.append(i)
            elif max_cd == f:
                max_cd_inds.append(i)
        if len(max_cd_inds) == len(fss):
            break
        for i in max_cd_inds:
            fs = fss[i]
            fs.pop(0)
            if len(fs) == 0:
                empty_fs_list += 1

    return max_cd

def normalize_ns(ns):
    """ Normlize a list of integers according to their gcd
    """
    gcd = gcd_brute(ns)
    return list(map(lambda v: v//gcd, ns))
def curdir():
    return os.path.abspath(os.curdir)

def load_cases_legacy(root):
    dir_stack = []
    cases = []
    dir_stack.append(curdir())
    os.chdir(root)
    for p in glob.iglob("**/profilecase_*", recursive=True):
        case_dir = abspath(p)
        dir_stack.append(curdir())
        os.chdir(case_dir)

        # Parse a case
        try:
            with open("profile.meta/meta.json", "r") as fh:
                case_entry = json.load(fh)
        except Exception as e:
            print("ERROR: {} misses meta.json".format(case_dir))
            os.chdir(dir_stack.pop())
            raise e

        timeout_init = False
        timeout_solve = False
        del case_entry["CommandTimeout"]
        del case_entry["Comments"]
        try:
            with open("profile.json", "r") as fh:
                case_profile = json.load(fh)
                case_entry["Machine"] = case_profile["info"]["machine"]
                del case_profile
        except Exception as e:
            timeout_init = os.path.exists("TIMEOUT_INIT")
            timeout_solve = os.path.exists("TIMEOUT_SOLVE")
            if not timeout_init and not timeout_solve:
                print("WARNING: {} misses profile.json".format(case_dir))
                os.chdir(dir_stack.pop())
                continue

        #print(case_dir)
        case_entry["CasePath"] = case_dir
        case_entry["TimeoutInit"] = timeout_init
        case_entry["TimeoutSolve"] = timeout_solve
        cases.append(case_entry)
        #print(case_entry)
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
    cases = load_cases_legacy(args.root_path)
    cases.to_csv(args.output_name)
