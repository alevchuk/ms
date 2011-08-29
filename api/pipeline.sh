set -e
set -u

set +u
if [ "$ALWAYS_RUN" == "" ]; then
  ALWAYS_RUN=false # Default is false, skip compleated steps whenever possible
fi
set -u

TIMEFORMAT="%0lR
" # When timing, print no-fraction, long-format, elapsed time, and a linebreak

function pipeline {

  current_prog=$0
  parent_prefix=$1  # Location(s) for the upstream time-stamp files
  prefix=$2         # Location for time-stamp files *-done and *-inprogress
  caption=$3
  payload=$4        # What needs to be done? (a callback function)

  echo $caption >&2

  in_progress=$prefix-inprogress-or-failed
  done_time=$prefix-done
  parent_time=$parent_prefix-done

  if ! $ALWAYS_RUN && [ ! -f $parent_prefix ] && [ ! -d $parent_prefix ]; then
    echo "ERROR: Upstream file or directory does not exist: $parent_prefix" >&2
    exit 1
  fi

  if ! $ALWAYS_RUN && [ ! -f $parent_time ]; then
    echo "ERROR: Upstream done file does not exist: $parent_time" >&2
    exit 1
  fi

  # Any one of these will trigger running the payload:
  #  * ALWAYS_RUN flag
  #  * Missig done file
  #  * My code file (current_prog) is newer than my done file
  #  * Parrent done file newer than my done file
  if $ALWAYS_RUN || [ ! -f $done_time ] || [ $current_prog -nt $done_time ] ||
    [ $parent_time -nt $done_time ]
  then
    time (
      mv $done_time $in_progress 2> /dev/null && true
      set +e
      touch $in_progress ||
        (echo "ERROR: Parent directory  must exist to track progress.";
        exit 1) || exit 1
      set -e
  
      $payload
  
      mv $in_progress $done_time
      touch $done_time
    )
  else
    printf "Already done.\n\n" >&2
  fi

}
