set -e
set -u

ALWAYS_RUN=false # Default should be true for portability
TIMEFORMAT="%0lR
" # When timing, print no-fraction, long-format, elapsed time, and a linebreak

function pipeline {

  parent_prefix=$1  # Location(s) for the upstream time-stamp files
  prefix=$2         # Location for time-stamp files *-done and *-inprogress
  caption=$3
  payload=$4        # What needs to be done? (a callback function)

  echo $caption 2> /dev/stderr

  in_progress=$prefix-inprogress-or-failed
  done_time=$prefix-done
  parent_time=$parent_prefix-done

  if ! $ALWAYS_RUN && [ ! -f $parent_time ]; then
    echo "ERROR: Parent done file does not exist: $parent_time" > /dev/stderr
    exit 1
  fi

  if $ALWAYS_RUN || [ ! -f $done_time ] || [ $parent_time -nt $done_time ]
  then
    time (
      mv $done_time $in_progress 2> /dev/null && true
      touch $in_progress
  
      $payload
  
      mv $in_progress $done_time
      touch $done_time
    )
  else
    printf "Already done.\n\n" 2> /dev/stderr
  fi

}
