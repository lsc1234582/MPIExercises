#!/bin/bash

TIME_OUT=30
SOLVER_TOLERANCE=0.00000001
X_MIN=-1.0
X_MAX=1.0
Y_MIN=-1.0
Y_MAX=1.0
# Test cases
declare -a FUNCTIONS=(0 1 2 3)
declare -a NUM_ROWS=(2 128)
declare -a NUM_COLS=(2 128)
declare -a NUM_PATCH_X=(1 2 3)
declare -a NUM_PATCH_Y=(1 2 3)

# Test stats
NUM_PASSES=0
NUM_TEST_CASES=0
SUMMARY_STR=""

# Get the source directory of the current script
# Source:
# https://stackoverflow.com/questions/59895/get-the-source-directory-of-a-bash-script-from-within-the-script-itself?page=1&tab=votes#tab-top
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# 1 - function selection
# 2 - num rows
# 3 - num cols
function test_serial () {
  test_case="SERIAL_FUNC${1}_NROW${2}_NCOL${3}"
  mkdir ${test_case}
  cd ${test_case}
  cat <<-SOMEMARK > parameters.in
XMin: ${X_MIN}
XMax: ${X_MAX}
YMin: ${Y_MIN}
YMax: ${Y_MAX}
NRow: ${2}
NCol: ${3}
Tolerance: ${SOLVER_TOLERANCE}
SOMEMARK
  echo "Running test case: ${test_case}"
  timeout ${TIME_OUT} ${DIR}/../Init parameters.in $1 2>&1 | tee logs
  timeout ${TIME_OUT} ${DIR}/../Solve parameters.in 2>&1 | tee -a logs
  python ${DIR}/test_checker.py -a solution.dat -b laplace.dat
  if [[ $? -eq 0 ]]; then
    echo "Test case: ${test_case} PASS"
    (( NUM_PASSES++ ))
    echo "PASS" > result.txt
    SUMMARY_STR+="${test_case}: PASS"$'\n'
  else
    echo "Test case: ${test_case} FAIL"
    echo "FAIL" > result.txt
    SUMMARY_STR+="${test_case}: FAIL"$'\n'
  fi
  (( NUM_TEST_CASES++ ))
  echo "========================================"
  cd ..
}

# 1 - function selection
# 2 - num rows
# 3 - num cols
# 4 - num patch in x
function test_parallel_row () {
  test_case="PARROW_FUNC${1}_NROW${2}_NCOL${3}_NUMPX${4}"
  mkdir ${test_case}
  cd ${test_case}
  cat <<-SOMEMARK > parameters.in
XMin: ${X_MIN}
XMax: ${X_MAX}
YMin: ${Y_MIN}
YMax: ${Y_MAX}
NRow: ${2}
NCol: ${3}
Tolerance: ${SOLVER_TOLERANCE}
SOMEMARK
  echo "Running test case: ${test_case}"
  timeout ${TIME_OUT} mpirun -n $4 ${DIR}/../InitPar $4 1 parameters.in $1 2>&1 | tee logs
  ${DIR}/../Combine $4 1 initial.combined.dat initial 2>&1 | tee -a logs
  ${DIR}/../Combine $4 1 solution.combined.dat solution 2>&1 | tee -a logs
  timeout ${TIME_OUT} mpirun -n $4 --mca btl self,tcp ${DIR}/../SolveParRow parameters.in 2>&1 | tee -a logs
  ${DIR}/../Combine $4 1 laplace.combined.dat laplace 2>&1 | tee -a logs
  python ${DIR}/test_checker.py -a solution.combined.dat -b laplace.combined.dat
  if [[ $? -eq 0 ]]; then
    echo "Test case: ${test_case} PASS"
    (( NUM_PASSES++ ))
    echo "PASS" > result.txt
    SUMMARY_STR+="${test_case}: PASS"$'\n'
  else
    echo "Test case: ${test_case} FAIL"
    echo "FAIL" > result.txt
    SUMMARY_STR+="${test_case}: FAIL"$'\n'
  fi
  (( NUM_TEST_CASES++ ))
  echo "========================================"
  cd ..
}

# 1 - function selection
# 2 - num rows
# 3 - num cols
# 4 - num patch in x
# 5 - num patch in y
function test_parallel () {
  test_case="PAR_FUNC${1}_NROW${2}_NCOL${3}_NUMPX${4}_NUMPY${5}"
  mkdir ${test_case}
  cd ${test_case}
  cat <<-SOMEMARK > parameters.in
XMin: ${X_MIN}
XMax: ${X_MAX}
YMin: ${Y_MIN}
YMax: ${Y_MAX}
NRow: ${2}
NCol: ${3}
Tolerance: ${SOLVER_TOLERANCE}
SOMEMARK
  echo "Running test case: ${test_case}"
  (( num_proc = $4 * $5 ))
  timeout ${TIME_OUT} mpirun -n ${num_proc} ${DIR}/../InitPar $4 $5 parameters.in $1 2>&1 | tee logs
  ${DIR}/../Combine $4 $5 initial.combined.dat initial 2>&1 | tee -a logs
  ${DIR}/../Combine $4 $5 solution.combined.dat solution 2>&1 | tee -a logs
  timeout ${TIME_OUT} mpirun -n ${num_proc} --mca btl self,tcp ${DIR}/../SolvePar $4 $5 parameters.in 2>&1 | tee -a logs
  ${DIR}/../Combine $4 $5 laplace.combined.dat laplace 2>&1 | tee -a logs
  python ${DIR}/test_checker.py -a solution.combined.dat -b laplace.combined.dat
  if [[ $? -eq 0 ]]; then
    echo "Test case: ${test_case} PASS"
    (( NUM_PASSES++ ))
    echo "PASS" > result.txt
    SUMMARY_STR+="${test_case}: PASS"$'\n'
  else
    echo "Test case: ${test_case} FAIL"
    echo "FAIL" > result.txt
    SUMMARY_STR+="${test_case}: FAIL"$'\n'
  fi
  (( NUM_TEST_CASES++ ))
  echo "========================================"
  cd ..
}

cd ${DIR}
if [ -d "build" ]; then
  rm -r build
fi
mkdir build
cd build

if [ -d "test_small" ]; then
  rm -r test_small
fi
mkdir test_small
cd test_small

# Run all serial test cases
for func in "${FUNCTIONS[@]}"
do
  for num_row in "${NUM_ROWS[@]}"
  do
    for num_col in "${NUM_COLS[@]}"
    do
      test_serial ${func} ${num_row} ${num_col}
    done
  done
done

# Run all parallel row test cases
for func in "${FUNCTIONS[@]}"
do
  for num_row in "${NUM_ROWS[@]}"
  do
    for num_col in "${NUM_COLS[@]}"
    do
      for num_patch_x in "${NUM_PATCH_X[@]}"
      do
        if (( ${num_row} >= ${num_patch_x} ))
        then
          test_parallel_row ${func} ${num_row} ${num_col} ${num_patch_x}
        fi
      done
    done
  done
done

# Run all parallel test cases
for func in "${FUNCTIONS[@]}"
do
  for num_row in "${NUM_ROWS[@]}"
  do
    for num_col in "${NUM_COLS[@]}"
    do
      for num_patch_x in "${NUM_PATCH_X[@]}"
      do
        for num_patch_y in "${NUM_PATCH_Y[@]}"
        do
          if (( ${num_row} >= ${num_patch_x} && ${num_col} >= ${num_patch_y} ))
          then
            test_parallel ${func} ${num_row} ${num_col} ${num_patch_x} ${num_patch_y}
          fi
        done
      done
    done
  done
done

# Print summary
echo "===== Test Summary ====="
echo "${SUMMARY_STR}"
echo "${NUM_PASSES} / ${NUM_TEST_CASES} passed"
