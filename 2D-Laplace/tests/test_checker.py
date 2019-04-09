import argparse
import numpy as np
import re
#matplotlib.use('agg')

# Tolerance params for checking discrepencies
ATOL = 1e-03
RTOL = 1e-02

def read_grid(fn):
    """Read and parse a grid and its parameters from a file. """
    params = {}
    grid = []
    row_ind = 0
    col_ind = 0
    x = 0
    y = 0
    x_min = None
    x_max = None
    y_min = None
    y_max = None
    num_row = None
    num_col = None
    is_last_match_value = False
    # File IO check
    # Syntactic formats correctness
    # Number of clumns in each row same
    # X coord monotinicity
    # Y coord monotinicity
    # X coord uniform increments
    # Y coord uniform increments
    with open(fn, "r") as fh:
        line = fh.readline()
        while line != "":
            match_obj = re.match("([\+-]?\d+\.?\d+)\s+([\+-]?\d+\.?\d+)\s+([\+-]?\d+\.?\d+)\s+", line)
            if match_obj is not None:
                # Matched a column
                col_ind += 1
                is_last_match_value = True
                x, y, val = [t(s) for t, s in zip((float, float, float), match_obj.groups())]
                if x_min is not None:
                    assert(x >= x_min)
                else:
                    x_min = x
                if y_min is not None:
                    assert(y >= y_min)
                else:
                    y_min = y
                grid.append(val)
            else:
                match_obj = re.match("\s*", line)
                if match_obj is not None:
                    # Matched a row
                    is_last_match_value = False
                    if num_col is None or num_col == col_ind:
                        num_col = col_ind
                        row_ind += 1
                        col_ind = 0
                    else:
                        raise Exception("Uneven number of columns across different rows")
                else:
                    raise Exception("Syntax error")
            line = fh.readline()

        # If the last matched line has values, then it means that there is an implicit delimiter (EOF)
        if is_last_match_value:
            if num_col is None or num_col == col_ind:
                num_col = col_ind
                row_ind += 1
                col_ind = 0
            else:
                raise Exception("Uneven number of columns across different rows")
        x_max = x
        y_max = y
        num_row = row_ind

    params["x_min"] = x_min
    params["x_max"] = x_max
    params["y_min"] = y_min
    params["y_max"] = y_max
    params["num_row"] = num_row
    params["num_col"] = num_col

    grid = np.array(grid).reshape((params["num_row"], params["num_col"]))

    return params, grid

def check_same(fn1, fn2):
    params1, grid1 = read_grid(fn1)
    params2, grid2 = read_grid(fn2)
    if params1 != params2:
        print("TEST_CHECKER: Parameters differ")
        return False
    diff = np.abs(grid2 - grid1)
    if not np.allclose(diff, 0.0, rtol=RTOL, atol=ATOL):
        print("TEST_CHECKER: Critical differences: mean {}, std {}".format(np.mean(diff), np.std(diff)))
        return False
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Laplace solver test checker")
    parser.add_argument("-a", "--grid_a_fn", dest="grid_a_fn", metavar="GridAFileName", action="store", type=str,
            help=("Grid A file name"), required=True)
    parser.add_argument("-b", "--grid_b_fn", dest="grid_b_fn", metavar="GridBFileName", action="store", type=str,
            help=("Grid B file name"), required=True)
    args = parser.parse_args()
    test_pass = False
    try:
        test_pass = check_same(args.grid_a_fn, args.grid_b_fn)
    except Exception as e:
        print("TEST_CHECKER: Caught unexpected exception")
        print(e)
        exit(1)
    if test_pass:
        exit(0)
    else:
        exit(1)
