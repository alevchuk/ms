TIMEFORMAT="%0lR
" # When timing, print no-fraction, long-format, elapsed time, and a linebreak

function pipeline {

  prefix=$1       # Location for the time-stamp (*-done and *-inprogress) files
  parrent_done=$2 # Location(s) for the upstream time-stamp files

  in_progress="$prefix-inprogress-or-failed"
  done_time="$prefix-done"
  if $ALWAYS_RUN || [ ! -f $done_time ] || [ $parrent_done -nt $done_time ]
  then
    time (
      mv $done_time $in_progress 2> /dev/null && true
      touch $in_progress
  
      payload
  
      mv $in_progress $done_time
      touch $done_time
    )
  else
    echo "Already done." 2> /dev/stderr
  fi

}
