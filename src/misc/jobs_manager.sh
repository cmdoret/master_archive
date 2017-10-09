# This script contains general utilities to manage LSF jobs.
# Cyril Matthey-Doret
# 10.10.2017


function bgenerate {
  # This function wraps a task to send it as a bjob
  # >>> Need to think of an efficient way to do this <<<
}

function bmonitor {
  # This function limits hangs the script when the number of queued bjobs
  # containing a given pattern in their namereaches a limit.
  # First argument is the name pattern of jobs to monitor
  # Second argument is the maximum of joobs that can be queued at a time.
  while [ $(bjobs -w | grep '$1' | wc -l) -gt $2 ]
  do
    sleep 1;
  done
}
