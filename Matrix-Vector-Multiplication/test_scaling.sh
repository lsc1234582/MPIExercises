NUM_REP=10

# Args
# $1 - Number of processes
# $2 - Number of repetitions
# $3 - ARow
# $4 - ACol
function record {
  local rep_count;
  rep_count=0;
  record_str="$1,$3,$4,";
  while [ ${rep_count} -le $2 ]; do
    record_str+="$(TIMEFORMAT="%3R"; time (mpirun -n $1 --oversubscribe ./main $3 $4) 2>&1),";
    rep_count=$(( rep_count + 1 ));
  done
  printf "${record_str}\n";
}

function test_strong_scaling {
  AROW=8192
  ACOL=8192
  # Number of processes = 1 (Serial)
  record 1 ${NUM_REP} ${AROW} ${ACOL};
  # Number of processes = 2
  record 2 ${NUM_REP} ${AROW} ${ACOL};
  # Number of processes = 4
  record 4 ${NUM_REP} ${AROW} ${ACOL};
  # Number of processes = 8
  record 8 ${NUM_REP} ${AROW} ${ACOL};
}

function test_weak_scaling {
  # Number of processes = 1 (Serial)
  AROW=2048
  ACOL=2048
  record 1 ${NUM_REP} ${AROW} ${ACOL};
  # Number of processes = 2
  AROW=4096
  ACOL=2048
  record 2 ${NUM_REP} ${AROW} ${ACOL};
  # Number of processes = 4
  AROW=8192
  ACOL=2048
  record 4 ${NUM_REP} ${AROW} ${ACOL};
  # Number of processes = 8
  AROW=16384
  ACOL=2048
  record 8 ${NUM_REP} ${AROW} ${ACOL};
}

# Main
cmd=$1
shift
${cmd} "$@"
