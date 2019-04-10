#!/bin/bash

TIME_OUT=120
SOLVER_TOLERANCE=0.00000001
X_MIN=-1.0
X_MAX=1.0
Y_MIN=-1.0
Y_MAX=1.0

# Profile cases
declare -a FUNCTIONS=(0)
#declare -a NUM_ROWS=(64 128 256 512)
#declare -a NUM_COLS=(64 128 256 512)
#declare -a NUM_PATCH_X=(1 2 3)
#declare -a NUM_PATCH_Y=(1 2 3)
declare -a NUM_ROWS=(64 128)
declare -a NUM_COLS=(64 128)
declare -a NUM_PATCH_X=(1 2)
declare -a NUM_PATCH_Y=(1 2)
## Version
PROFILE_SERIAL=0
PROFILE_PARROW=0
PROFILE_PAR=1
BUILD_CMD="make DEBUG=1"
CLEAN_CMD="make clean"
SOLVER_LAUNCH_CMD="mpirun --mca btl self,tcp"
REP=3 # Number of repetitions for each profile case

## Code version and comment: automatically save git log and git branch output, and associate a comment with the code
    ## version
ADDITINAL_TAGS=""
COMMENTS=""

# Profile stats
SUMMARY_STR=""
NUM_PROFILE_CASES=0

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

# 1 - NRow
# 2 - NCol
function write_parameters () {
  cat <<-SOMEMARK > parameters.in
XMin: ${X_MIN}
XMax: ${X_MAX}
YMin: ${Y_MIN}
YMax: ${Y_MAX}
NRow: ${1}
NCol: ${2}
Tolerance: ${SOLVER_TOLERANCE}
SOMEMARK
}

function build () {
  pushd ${DIR}/..
  ${CLEAN_CMD}
  ${BUILD_CMD}
  popd
}

# 1 - Path to profile.meta
# 2 - StartDateTime
# 3 - EndDateTime
# 4 - SolverLaunchCommand
# 5 - Solver
# 6 - NRow
# 7 - NCol
# 8 - NPX
# 9 - NPY
function generate_meta_json () {
  git branch | grep \* | cut -d ' ' -f2 > $1/git.info
  git log --oneline >> $1/git.info
  git diff >> $1/git.info
  git_tag=$(git tag --points-at HEAD)
  if [ "${git_tag}" = "" ]; then
    echo "Error: Current commit must be tagged"
    exit 1
  fi
  cat <<-SOMEMARK > $1/meta.json
{
"StartDateTime": "$2",
"EndDateTime": "$3",
"CommandTimeout" : "${TIME_OUT}",
"Solver": "$5",
"NRow": "$6",
"NCol": "$7",
"NPX": "$8",
"NPY": "$9",
"SourceVersionTag": "${git_tag}",
"BuildCommand": "${BUILD_CMD}",
"SolverLaunchCommand": "$4",
"AdditionalTags": "${ADDITINAL_TAGS}",
"Comments": "${COMMENTS}"
}
SOMEMARK
}

# 1 - function selection
# 2 - num rows
# 3 - num cols
function profile_serial () {
  profile_case="SERIAL_FUNC${1}_NROW${2}_NCOL${3}_NPX${4}"
  date=$(date +"D%dM%mY%yH%HM%MS%S")
  profile_case_dir="profilecase_${profile_case}_${date}"
  mkdir ${profile_case_dir}
  cd ${profile_case_dir}
  start_date_time=$(date)
  mkdir profile.meta
  write_parameters $2 $3
  echo "Building profile case: ${profile_case}"
  build 2>&1 > profile.meta/compile.info
  echo "Running profile case: ${profile_case}"
  timeout ${TIME_OUT} ${DIR}/../Init parameters.in $1 2>&1 > logs
  timeout ${TIME_OUT} map --profile --log=map_profile_log --output=profile --nompi ${DIR}/../Solve parameters.in 2>&1 >> logs
  map --profile --export=profile.json profile.map
  echo "========================================"
  end_date_time=$(date)
  generate_meta_json "$(pwd)/profile.meta" "${start_date_time}" "${end_date_time}" "" "Serial" "$2" "$3" "" ""
  (( NUM_PROFILE_CASES++ ))
  cd ..
}

# 1 - function selection
# 2 - num rows
# 3 - num cols
# 4 - num patch in x
function profile_parallel_row () {
  profile_case="PARROW_FUNC${1}_NROW${2}_NCOL${3}_NPX${4}"
  date=$(date +"D%dM%mY%yH%HM%MS%S")
  profile_case_dir="profilecase_${profile_case}_${date}"
  mkdir ${profile_case_dir}
  cd ${profile_case_dir}
  start_date_time=$(date)
  mkdir profile.meta
  write_parameters $2 $3
  echo "Building profile case: ${profile_case}"
  build 2>&1 > profile.meta/compile.info
  echo "Running profile case: ${profile_case}"
  timeout ${TIME_OUT} mpirun -n $4 ${DIR}/../InitPar $4 1 parameters.in $1 2>&1 > logs
  ${DIR}/../Combine $4 1 initial.combined.dat initial 2>&1 >> logs
  ${DIR}/../Combine $4 1 solution.combined.dat solution 2>&1 >> logs
  timeout ${TIME_OUT} map --profile --log=map_profile_log --output=profile ${SOLVER_LAUNCH_CMD} -n $4 ${DIR}/../SolveParRow parameters.in 2>&1 >> logs
  ${DIR}/../Combine $4 1 laplace.combined.dat laplace 2>&1 >> logs
  map --profile --export=profile.json profile.map
  echo "========================================"
  end_date_time=$(date)
  generate_meta_json "$(pwd)/profile.meta" "${start_date_time}" "${end_date_time}" "${SOLVER_LAUNCH_CMD}" "ParRow" "$2" "$3" "$4" ""
  (( NUM_PROFILE_CASES++ ))
  cd ..
}

# 1 - function selection
# 2 - num rows
# 3 - num cols
# 4 - num patch in x
# 5 - num patch in y
function profile_parallel () {
  profile_case="PAR_FUNC${1}_NROW${2}_NCOL${3}_NPX${4}_NPY${5}"
  date=$(date +"D%dM%mY%yH%HM%MS%S")
  profile_case_dir="profilecase_${profile_case}_${date}"
  mkdir ${profile_case_dir}
  cd ${profile_case_dir}
  start_date_time=$(date)
  mkdir profile.meta
  write_parameters $2 $3
  echo "Building profile case: ${profile_case}"
  build 2>&1 > profile.meta/compile.info
  echo "Running profile case: ${profile_case}"
  (( num_proc = $4 * $5 ))
  timeout ${TIME_OUT} mpirun -n ${num_proc} ${DIR}/../InitPar $4 $5 parameters.in $1 2>&1 > logs
  ${DIR}/../Combine $4 $5 initial.combined.dat initial 2>&1 >> logs
  ${DIR}/../Combine $4 $5 solution.combined.dat solution 2>&1 >> logs
  timeout ${TIME_OUT} map --profile --log=map_profile_log --output=profile ${SOLVER_LAUNCH_CMD} -n ${num_proc} ${DIR}/../SolvePar $4 $5 parameters.in 2>&1 >> logs
  ${DIR}/../Combine $4 $5 laplace.combined.dat laplace 2>&1 >> logs
  map --profile --export=profile.json profile.map
  echo "========================================"
  end_date_time=$(date)
  generate_meta_json "$(pwd)/profile.meta" "${start_date_time}" "${end_date_time}" "${SOLVER_LAUNCH_CMD}" "Par" "$2" "$3" "$4" "$5"
  (( NUM_PROFILE_CASES++ ))
  cd ..
}

cd ${DIR}
git_tag=$(git tag --points-at HEAD)
date=$(date +"D%dM%mY%yH%HM%MS%S")
build_dir="profile_${git_tag}_${date}"
if [ -d ${build_dir} ]; then
  echo "Error: same build folder already exists"
  exit 1
fi
mkdir ${build_dir}
cd ${build_dir}

if [[ ${PROFILE_SERIAL} -eq 1 ]]; then
  # Run all serial profile cases
  for func in "${FUNCTIONS[@]}"
  do
    for num_row in "${NUM_ROWS[@]}"
    do
      for num_col in "${NUM_COLS[@]}"
      do
        for (( i=0; i<${REP}; i++));
        do
          profile_serial ${func} ${num_row} ${num_col}
        done
      done
    done
  done
fi

# Run all parallel row profile cases
if (( ${PROFILE_PARROW} == 1 )); then
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
            for (( i=0; i<${REP}; i++));
            do
              profile_parallel_row ${func} ${num_row} ${num_col} ${num_patch_x}
            done
          fi
        done
      done
    done
  done
fi

# Run all parallel profile cases
if (( ${PROFILE_PAR} == 1 )); then
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
              for (( i=0; i<${REP}; i++));
              do
                profile_parallel ${func} ${num_row} ${num_col} ${num_patch_x} ${num_patch_y}
              done
            fi
          done
        done
      done
    done
  done
fi

# Print summary
echo "===== Profile Summary ====="
echo "${SUMMARY_STR}"
echo "${NUM_PROFILE_CASES} profile cases run"
