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
import subprocess as spr


# Solver category order mapping: Axis type -> key func
solver_map = {k:v + 1 for v, k in enumerate(["Serial", "ParRow", "Par"])}
func_map = {k:v + 1 for v, k in enumerate(["0", "1", "2", "3"])}

int_axes = ["NRow", "NCol", "NPX", "NPY", "SampleSize"]
float_axes = ["ComputeTime", "IOTime", "MPITime", "SolverTolerance", "XMax", "XMin", "YMax", "YMin"]
str_axes = ["AdditionalTags", "BuildCommand", "Func", "Machine", "Solver", "SolverLaunchCommand", "SourceVersionTag"]
all_all_index_axes = ["Machine", "SourceVersionTag", "BuildCommand",
                  "Func", "Size", "NRow", "NCol", "SolverTolerance", "XMax", "XMin", "YMax", "YMin",
                  "Solver", "Proc", "NPX", "NPY",
                 "AdditionalTags",
                 "StartDateTime", "EndDateTime", "TimeoutInit", "TimeoutSolve"]

all_index_axes = ["Machine", "SourceVersionTag", "BuildCommand",
                  "Func", "Size", "NRow", "NCol", "SolverTolerance", "XMax", "XMin", "YMax", "YMin",
                "AdditionalTags",
                "TimeoutInit", "TimeoutSolve",
                "Solver", "Proc", "NPX", "NPY"]

all_index_axes_final = ["Machine", "SourceVersionTag", "BuildCommand",
                  "Func", "Size", "NRow", "NCol", "SolverTolerance", "XMax", "XMin", "YMax", "YMin",
                 "AdditionalTags",
                "TimeoutInit", "TimeoutSolve",
                "Solver", "Proc", "NPX", "NPY",
                       "RepID"]

composite_axes = ["Size:NRow^NCol", "Proc:NPX^NPY"]
categorical_axes = ["AdditionalTags", "BuildCommand", "Func", "Machine", "Solver", "SolverLaunchCommand", "SourceVersionTag",
        "TimeoutInit", "TimeoutSolve"]


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

# profile stats

################################ Utilities ##############################

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

def is_floating(t):
    return t in [float, np.float32, np.float64]

def get_index_levels(df, index_name, axis=0):
    """ Get a list of levels for a index in a potentially multi-index df
    """
    if axis == 0:
        index = df.index
    elif axis == 1:
        index = df.columns
    else:
        assert(False)

    # The numeric index of the index associated with index_name
    index_idx = index.names.index(index_name)
    return index.levels[index_idx]

################################### Case processing ######################################

def add_composite_axes(cases, composite_axes):
    # Composite axes : composite axis name -> func
    cases = cases.copy()
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
        if is_integral(cases.loc[:, axes[0]].dtype) and is_integral(cases.loc[:, axes[1]].dtype):
            c_axes[c_ax_name] = (lambda subaxes: (lambda x: x[subaxes[0]] * x[subaxes[1]]))(axes)
        else:
            c_axes[c_ax_name] = (lambda subaxes: (lambda x: str(x[subaxes[0]]) + "^" + str(x[subaxes[1]])))(axes)

    return cases.assign(**c_axes)

def remove_useless_axes(cases):
    return cases.drop(columns=["SolverLaunchCommand", "EndDateTime", "StartDateTime"])


def add_repid(cases, all_index_axes):
    # Add RepID
    for ind, new_df in cases.groupby(all_index_axes):
        rep_id = 1
        for i, row in new_df.iterrows():
            cases.at[i, "RepID"] = int(rep_id)
            rep_id += 1
    return cases.astype({"RepID":"int"})

def index_axes(cases, all_index_axes):
    cases.set_index(all_index_axes, inplace=True)
    cases.sort_index(inplace=True)
    return cases

def create_unstack_cases(cases):
    cases_unstacked = cases.unstack(["Solver", "Proc", "NPX", "NPY"]).reorder_levels(["Solver", "Proc", "NPX", "NPY", 0], axis=1).sort_index(axis=1)
    cases_unstacked.columns.names=["Solver", "Proc", "NPX", "NPY", "Metrics"]
    return cases_unstacked

#################################### Case enhancement ######################################
def calculate_par_solver_speedup(cases_unstacked):
    # Calculate parallel speed up (per solver) for Par solvers
    num_npx = len(get_index_levels(cases_unstacked, "NPX", 1))
    num_npy = len(get_index_levels(cases_unstacked, "NPY", 1))
    serial_runtime = cases_unstacked.xs(["Par", "RunTime", 1, 1, 1], level=["Solver", "Metrics", "Proc", "NPX", "NPY"], axis=1)
    ind = cases_unstacked.xs(["Par", "RunTime"], level=["Solver", "Metrics"], axis=1, drop_level=False).columns
    serial_runtime_dup = pd.concat([serial_runtime for _ in range(num_npx * num_npy)], axis=1)
    serial_runtime_dup.columns = ind
    cases_speedup_pers = serial_runtime_dup / cases_unstacked.xs(["Par", "RunTime"], level=["Solver", "Metrics"], axis=1, drop_level=False)

    cases_speedup_pers = cases_speedup_pers.rename({"RunTime": "ParSpeedupPerSolver"}, axis=1)
    return cases_speedup_pers

def calculate_par_solver_efficiency(cases_speedup_pers):
    # Calculate efficiency (per solver) for Par solvers
    cases_unstacked_all_metrics_copy = cases_speedup_pers.xs("ParSpeedupPerSolver", level="Metrics", axis=1, drop_level=False).copy()
    ind = cases_unstacked_all_metrics_copy.index
    col = cases_unstacked_all_metrics_copy.columns
    proc_ind = cases_unstacked_all_metrics_copy.columns.names.index("Proc")
    proc_codes = cases_unstacked_all_metrics_copy.columns.codes[proc_ind]
    procs = [cases_unstacked_all_metrics_copy.columns.levels[proc_ind][c] for c in proc_codes]
    procs = pd.DataFrame([procs for _ in range(len(cases_unstacked_all_metrics_copy))], index=ind, columns=col)
    cases_efficiencies_pers = (cases_unstacked_all_metrics_copy / procs).rename({"ParSpeedupPerSolver": "ParEfficiencyPerSolver"}, axis=1)
    return cases_efficiencies_pers

def calculate_runtimemultiplier(cases_unstacked, all_index_axes_final, solver):
    cases_size_unstacked = cases_unstacked.stack(["Solver", "Proc", "NPX", "NPY"]).reorder_levels(all_index_axes_final).sort_index(axis=1).reset_index().copy()
    cases_size_unstacked.set_index(all_index_axes_final, inplace=True)
    cases_size_unstacked.sort_index(inplace=True)
    cases_size_unstacked = cases_size_unstacked.unstack(["Solver", "Size", "NRow", "NCol"]).reorder_levels(["Solver", "Size", "NRow", "NCol", "Metrics"], axis=1).sort_index(axis=1)
    cases_size_unstacked.columns.names=["Solver", "Size", "NRow", "NCol", "Metrics"]

    size_levels = get_index_levels(cases_size_unstacked, "Size", 1)
    basesize = size_levels[0]

    nrow_levels = get_index_levels(cases_size_unstacked, "NRow", 1)
    basenrow = nrow_levels[0]

    ncol_levels = get_index_levels(cases_size_unstacked, "NCol", 1)
    basencol = ncol_levels[0]

    cases_runtimemult = []
    #print(len(size_levels))
    basesize_runtime = cases_size_unstacked.xs([solver, "RunTime", basesize, basenrow, basencol], level=["Solver", "Metrics", "Size", "NRow", "NCol"], axis=1)
    ind = cases_size_unstacked.xs([solver, "RunTime"], level=["Solver", "Metrics"], axis=1, drop_level=False).columns
    basesize_runtime_dup = pd.concat([basesize_runtime for _ in range(len(nrow_levels) * len(ncol_levels))], axis=1)
    basesize_runtime_dup.columns = ind
    cases_runtimemult.append(cases_size_unstacked.xs([solver, "RunTime"], level=["Solver", "Metrics"], axis=1, drop_level=False) / basesize_runtime_dup)

    #cases_size_unstacked.xs(["RunTime"], level=["Metrics"], axis=1, drop_level=False)
    cases_runtimemult = pd.concat(cases_runtimemult, axis=1).rename({"RunTime": "RunTimeMultiplier"}, axis=1)

    cases_runtimemult = cases_runtimemult.stack(["Solver", "Size", "NRow", "NCol"]).reorder_levels(all_index_axes_final).sort_index(axis=1).reset_index().copy()
    cases_runtimemult.set_index(all_index_axes_final, inplace=True)
    cases_runtimemult.sort_index(inplace=True)
    cases_runtimemult = cases_runtimemult.unstack(["Solver", "Proc", "NPX", "NPY"]).reorder_levels(["Solver", "Proc", "NPX", "NPY", "Metrics"], axis=1).sort_index(axis=1)
    cases_runtimemult.columns.names=["Solver", "Proc", "NPX", "NPY", "Metrics"]
    return cases_runtimemult


def enhance_serial_cases(cases_unstacked, all_index_axes_final):
    cases_runtimemult = calculate_runtimemultiplier(cases_unstacked, all_index_axes_final, "Serial")
    cases_unstacked_all_metrics = pd.concat([cases_unstacked,
                                            cases_runtimemult], axis=1, sort=True)
    return cases_unstacked_all_metrics


def enhance_par_cases(cases_unstacked, all_index_axes_final):
    cases_speedup_pers = calculate_par_solver_speedup(cases_unstacked)
    cases_efficiencies_pers = calculate_par_solver_efficiency(cases_speedup_pers)
    cases_runtimemult = calculate_runtimemultiplier(cases_unstacked, all_index_axes_final, "Par")
    cases_unstacked_all_metrics = pd.concat([cases_unstacked,
                                            cases_speedup_pers,
                                            cases_efficiencies_pers,
                                            cases_runtimemult], axis=1, sort=True)
    return cases_unstacked_all_metrics


def select_cases(cases, sel, metrics, legacy_load=False, copy_cases=False):
    #group_by_axes = axes[0].split("^") if len(axes) == 1 and "^" in axes[0] else axes
    if copy_cases:
        plot_df = cases.copy()
    else:
        plot_df = cases

    if len(sel) != 0:
        plot_df = plot_df.query(sel)

    if legacy_load:
        for i, row in plot_df.iterrows():
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

    if "Speedup" in metrics:
        plot_df = calculate_speedup(plot_df)
    if "RunTimeNorm" in metrics:
        plot_df = calculate_runtimenorm(plot_df)

    return plot_df

def calculate_speedup(df):
    max_runtime = np.max(df.loc[:, "RunTime"])
    df = df.assign(Speedup=lambda x: max_runtime / x["RunTime"])
    return df

def calculate_runtimenorm(df):
    min_runtime = np.min(df["RunTime"])
    df = df.assign(RunTimeNorm=lambda x: x["RunTime"] / min_runtime)
    return df

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
    """ Yield all profile cases, one by one, from all profiles (directory and .tar.gz) at root folder.
    If it's a .tar.gz archive, uncompress to uncompress_folder. This folder will be deleted if delete_uncompress_folder
    is True.
    Furthermore, if uncompress_folder already exists an exception will be raised.
    """

    dir_stack = []
    dir_stack.append(curdir())
    os.chdir(root)
    if os.path.exists(uncompress_folder):
        raise Exception("Uncompress folder {} already exists".format(uncompress_folder))
    print("Making temporary uncompress folder {}...".format(uncompress_folder))
    os.makedirs(uncompress_folder)

    def is_tar_gz_file(f):
        return os.path.isfile(f) and tarfile.is_tarfile(f) and "gz" in os.path.splitext(f)[1]

    for p in glob.iglob("**/profile_*", recursive=True):
        print("INFO: profile folder/file {} detected".format(p))
        if p.split(os.path.sep)[0] == uncompress_folder:
            # Make sure do not re-glob uncompressed folders
            print("Profile {} already exists in uncompress_folder {}. Skip uncompressing".format(p, uncompress_folder))
            continue
        elif os.path.isdir(p):
            print("Profile {} is a directory".format(p))
            profile_dir = abspath(p)
        elif is_tar_gz_file(p):
            print("Profile {} is a compressed tar ball".format(p))
            p_without_ext = os.path.splitext(p)[0]
            p_without_ext = os.path.splitext(p_without_ext)[0]
            tmp_p_dir = os.path.join(uncompress_folder, p_without_ext)
            os.makedirs(tmp_p_dir)
            profile_dir = abspath(tmp_p_dir)
            with tarfile.open(p) as t:
                if uncompress_member_func is None:
                    t.extractall(path=profile_dir)
                else:
                    t.extractall(members=uncompress_member_func(t), path=profile_dir)
        else:
            print("Profile {} is not recognised. Skipping...".format(p))
            continue

        dir_stack.append(curdir())
        os.chdir(profile_dir)
        yield from glob.iglob("**/profilecase_*", recursive=True)
        os.chdir(dir_stack.pop())

        if delete_uncompress_folder and is_tar_gz_file(p):
            print("Deleting temporary profile folder {}...".format(profile_dir))
            sh.rmtree(profile_dir, ignore_errors=True)
    if delete_uncompress_folder:
        print("Deleting temporary uncompress folder {}...".format(uncompress_folder))
        sh.rmtree(uncompress_folder, ignore_errors=True)

def handle_error_string(fn, err_pattern_str, rpl_str):
    err_pattern = re.compile(err_pattern_str)
    with open(fn, "r") as fh:
        lines = fh.read()
    match = err_pattern.search(lines)
    if match is None:
        return False
    correct_str = err_pattern.sub(rpl_str, lines)
    with open(fn, "w") as fh:
        fh.write(correct_str)
    return True

def load_cases(root, metrics):
    """Discover all profile cases under a root dir. The profile cases can be in a folder or a tar.gz compressed archive
    Return a DataFrame object containing all the params combination and a path to each case.
    """
    dir_stack = []
    cases = []
    timed_out_cases = []

    def profile_and_meta_json_file(members):
        for tarinfo in members:
            if "profile_" in tarinfo.name and\
                    "profilecase_" in tarinfo.name and\
                    ("profile.meta/meta.json" in tarinfo.name or\
                     "profile.json" in tarinfo.name):
                yield tarinfo

    for case_dir in generate_profile_cases(root, ".profiletmp", True):
        dir_stack.append(curdir())
        os.chdir(case_dir)
        # Parse a case
        try:
            with open("profile.meta/meta.json", "r") as fh:
                case_entry = json.load(fh)
        except Exception as e:
            missing_npx = handle_error_string("profile.meta/meta.json",
                    r'"NPX":\s*,',
                    r'"NPX": 1,')
            missing_npy = handle_error_string("profile.meta/meta.json",
                    r'"NPY":\s*,',
                    r'"NPY": 1,')
            if missing_npx or missing_npy:
                with open("profile.meta/meta.json", "r") as fh:
                    case_entry = json.load(fh)
            else:
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
            if timeout_init or timeout_solve:
                print("WARNING: {} timed out, skipping".format(case_dir))
                timed_out_cases.append(case_dir)
                os.chdir(dir_stack.pop())
                continue
            if os.path.exists("profile.map"):
                print("INFO: {} misses profile.json but found profile.map. Generating profile.json".format(case_dir))
                #spr.run(["module", "load", "arm-hpc-tools/arm-forge"], shell=True, executable="/bin/bash")
                spr.run(["map", "--profile", "--export=profile.json", "profile.map"])
                with open("profile.json", "r") as fh:
                    profile_json = json.load(fh)
            else:
                print("WARNING: {} misses profile.json and profile.map".format(case_dir))
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
    print(timed_out_cases)
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

def post_load_case(cases, legacy=False):
    """ Post-processing after loading cases
    TODO: Get rid of errorneous and inrecoverable cases.
    """
    cases = cases.copy()
    # Fill in missing values (NPX, NPY)
    # if legacy:
    #     cases.loc[(cases.loc[:, "NPX"]) == "", "NPX"] = 1
    #     cases.loc[(cases.loc[:, "NPY"]) == "", "NPY"] = 1
    #     cases.loc[:, "machine"] = cases.loc[:, "machine"].fillna("e121008-lin")
    #     cases = cases.assign(SolverTolerance=0.00000001,
    #             XMin=-1.0,
    #             XMax=1.0,
    #             YMin=-1.0,
    #             YMax=1.0)
    # # Transform type
    # if legacy:
    #     cases.loc[:, ["NCol", "NRow", "NPX", "NPY"]] = cases.loc[:, ["NCol", "NRow", "NPX", "NPY"]].applymap(int)
    # # Make categories
    return cases


def correct_dtypes(cases):
    """ Pre-process: correct data types of the cases
    """
    cases = cases.copy()
    #cases = cases.drop(["Unnamed: 0"], axis=1)

    cases.loc[:, int_axes] = cases.loc[:, int_axes].fillna(0)
    cases.loc[:, float_axes] = cases.loc[:, float_axes].fillna(0.0)
    cases.loc[:, str_axes] = cases.loc[:, str_axes].fillna("")

    cases.loc[:, int_axes] = cases.loc[:, int_axes].astype("int")
    cases.loc[:, float_axes] = cases.loc[:, float_axes].astype("float")
    cases.loc[:, str_axes] = cases.loc[:, str_axes].astype("str")

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
    timed_out_cases = []
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
            if timeout_init or timeout_solve:
                timed_out_cases.append(case_dir)
            if not timeout_init and not timeout_solve:
                if os.path.exists("profile.map"):
                    print("INFO: {} misses profile.json but found profile.map. Generating profile.json".format(case_dir))
                    #spr.run(["module", "load", "arm-hpc-tools/arm-forge"], shell=True, executable="/bin/bash")
                    #spr.run(["map", "--profile", "--export=profile.json", "profile.map"], shell=True,
                    #        executable="/bin/bash")
                    #with open("profile.json", "r") as fh:
                    #    case_profile = json.load(fh)
                    #    case_entry["Machine"] = case_profile["info"]["machine"]
                    #    del case_profile
                    os.chdir(dir_stack.pop())
                    continue
                else:
                    print("WARNING: {} misses profile.json and profile.map".format(case_dir))
                    os.chdir(dir_stack.pop())
                    continue

        #print(case_dir)
        case_entry["CasePath"] = case_dir
        case_entry["TimeoutInit"] = timeout_init
        case_entry["TimeoutSolve"] = timeout_solve
        cases.append(case_entry)
        #print(case_entry)
        os.chdir(dir_stack.pop())

    print(timed_out_cases)
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
