#!/usr/bin/env ruby

msa_list = Hash.new []
while line = STDIN.gets
  msa_name, seq_id, seq = line.split("\t")
  seq.gsub!("-", "")
  msa_list[msa_name] = msa_list[msa_name] + [seq]
end

msa_list.each_pair do |msa, seq_list|
  numseq = seq_list.size
  numres = seq_list.inject(0){|accum,y| accum + y.size}
  avrseqlen = (numres / numseq).round # Try to round consitenty with
                                      # 004-reshape-cdd-scripts/filter-cdd.R
  puts "#{msa}\t#{numseq}\t#{numres}\t#{avrseqlen}"
end
