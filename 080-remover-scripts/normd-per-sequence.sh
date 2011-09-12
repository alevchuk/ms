#!/bin/bash

# UNNAMED is a classifier of outlier sequences in multiple alignments
# Copyright (C) 2011  Aleksandr Levchuk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


set -e
set -u

if [ $# != 2 ]; then
  echo "Usage: $0 <opt-dir> <remover-dir>" >&2
  echo "" >&2
  echo "  <opt-dir>        where dependeny softwre is installed" >&2
  echo "  <remover-dir>    where everithing is happening" >&2
  echo "" >&2
  exit 1
fi

# Arguments
OPT_DIR=$1
REMOVER_DIR=$2

# This is the output file!
out_dx1_mod_file=$REMOVER_DIR/normd-dx1-mod-seq-scores.txt

RAM_DIR=/dev/shm/`basename $REMOVER_DIR`

NORMD="$OPT_DIR/norMD1_3/normd"
SCORING_MATRIX="$OPT_DIR/norMD1_3/gon250.bla"
NORMD_OPTS="$SCORING_MATRIX 0 0 0 -v"
MAFFT_PATH="$OPT_DIR/mafft-6.857-with-extensions/bin"
CONVERT_FASTA_TO_TAB="$OPT_DIR/misc-bioinfo-scripts/convert-fasta-to-tab"

seq_group_fasta=$REMOVER_DIR/msa-fasta-tab-nogaps-fasta.txt
seq_group_tab=$REMOVER_DIR/msa-fasta-tab-nogaps.txt
num_seq=$(wc -l $seq_group_tab| awk '{print $1}')

mkdir -p $RAM_DIR 2> /dev/null && true

# Align
msa_fasta=$REMOVER_DIR/msa-mafft-fasta.txt
$MAFFT_PATH/mafft  --reorder --amino  --quiet $seq_group_fasta > $msa_fasta
msa_tab=$msa_fasta-tab.txt
cat $msa_fasta | $CONVERT_FASTA_TO_TAB > $msa_tab


# Get the whole MSA score
set +e
score_raw=$($NORMD $msa_fasta $NORMD_OPTS)
set -e
(echo "$score_raw" | grep "ERROR") && (echo "$score_raw" >&2 exit 1 )
whole_normd=$(echo "$score_raw" | egrep "^norMD" | awk '{print $2}')
whole_lqrid=$(echo "$score_raw" | egrep "^LQR" | awk '{print $2}')
whole_md=$(echo "$score_raw" | egrep "^MD" | awk '{print $2}')
whole_maxmd=$(echo "$score_raw" | egrep "^maxMD" | awk '{print $2}')
WHOLE_NORMD_MOD=$(ruby -e "puts $whole_md / $whole_maxmd")


# Remove sequences and score
ram_file1=$RAM_DIR/normd-dx1-mod-seq-scores.txt
cat /dev/null > $ram_file1
for seq_id in `seq 1 $num_seq`; do
  ram_file2=$RAM_DIR/seq$seq_id-removed
  seq_name=$(cat $msa_tab | awk "NR==$seq_id" | awk '{print $1}')

  # Remove a sequence
  cat $msa_tab | awk "NR!=$seq_id" > $ram_file2

  # Convert tab to FASTA
  cat $ram_file2 | awk '{print ">" $1 "\n" $2}' > $ram_file2-fasta

  # Run NorMD 
  set +e
  score_raw=$($NORMD $ram_file2-fasta $NORMD_OPTS)
  set -e
  (echo "$score_raw" | grep "ERROR") && (echo "$score_raw" >&2 exit 1)

  normd=$(echo "$score_raw" | egrep "^norMD" | awk '{print $2}')
  lqrid=$(echo "$score_raw" | egrep "^LQR" | awk '{print $2}')
  md=$(echo "$score_raw" | egrep "^MD" | awk '{print $2}')
  maxmd=$(echo "$score_raw" | egrep "^maxMD" | awk '{print $2}')
  # By the way: NORMD = MD / (MaxMD * LQRID)

  normd_dx1_mod=$(ruby -e "puts $WHOLE_NORMD_MOD - ($md / $maxmd)")

  #### This as well as the absolute $normd results in random guessing
  #modified_normd=$(ruby -e "puts $md / $maxmd")

  # Output
  echo -e "$seq_name\t$normd_dx1_mod" >> $ram_file1

  rm $ram_file2
done

# Can't use sort -n because the file is in scientific notation
cat $ram_file1 | ruby -e 'STDIN.readlines.collect{|x|
  v1,v2=x.split
  [v1, v2.to_f]
}.sort{|i1,i2|
  i1[1] <=> i2[1]}.each{|x| puts "#{x[0]}\t#{x[1]}"}' > $ram_file1-srt
mv $ram_file1-srt   $out_dx1_mod_file

rmdir $RAM_DIR 2> /dev/null && true

exit 0
