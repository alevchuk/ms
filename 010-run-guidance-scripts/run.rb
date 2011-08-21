#!/usr/bin/env ruby

# Limit to avoid waiting forever
MAX_FASTA_SIZE = 200_000 # number of non-header characters

RUNDIR = "011-run1-4-mafft-v6.712"

Dir.mkdir(RUNDIR) unless File.directory?(RUNDIR)

fasta_files = []
Dir.glob("cw/20*/*.fasta") do |fname|

  s = `cat #{fname} | grep -v "^>" | wc -c`.to_i

  if s > MAX_FASTA_SIZE
    STDERR.puts "Skipping because fasta size > #{MAX_FASTA_SIZE}:" +
      " #{fname} (#{s} characters)"
    next
  end

  fasta_files.push fname
end


fasta_files.each do |fpath|
  if fpath =~ /\/([^\/]*)\/([^\/]*)\.fasta/
    simple_name = "#{$1}.#{$2}"
  
    # --msaParam --thread 8 
    opts = "--msaProgram MAFFT --seqType aa"
    cmd = "guidance.pl #{opts} " + 
      "--seqFile ../#{fpath}  --outDir #{simple_name}.guidance"

    File.open("./#{RUNDIR}/#{simple_name}.sh", 'w') do |f|
      f.puts "export PATH=~/opt/mafft-6.712-with-extensions/bin:$PATH"
      f.puts "(which mafft guidance.pl; time #{cmd}) &> " +
        "./#{simple_name}.sh.out"
    end

    qsub_cmd = "qsub -d #{Dir.pwd}/#{RUNDIR} ./#{simple_name}.sh"

    puts qsub_cmd
  else
    raise
  end

end