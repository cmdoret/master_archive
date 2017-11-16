# This script contains general utilities to manage LSF jobs.
# Cyril Matthey-Doret
# 10.10.2017



function bmonitor {
  # This function limits hangs the script when the number of queued bjobs
  # containing a given pattern in their namereaches a limit.
  # First argument is the name pattern of jobs to monitor
  # Second argument is the maximum of joobs that can be queued at a time.
  while [ $(bjobs -w | grep "$1" | wc -l) -gt $2 ]
  do
    sleep 1;
  done
}

function prettyload {
  # This function displays a small animation and shows real-time
  # progress of a task.
  case $(expr $1 % 4) in
    0) sym="|";;
    1) sym="/";;
    2) sym="-";;
    3) sym="\\";;
  esac
  echo -ne "  $sym $1/$2 \r"
}
