%%%
%%% createrPBSfile_Hoffman.m
%%%
%%% Writes out a script that runs the MITgcm. Takes as arguments the name
%%% of the experiment, the number of nodes that will be used, and the wall
%%% time in hours (after which the computation will be cut off regardless 
%%% of whether it is complete).
%%%
function createPBSfile_Hoffman (dirname,expname,nodes)

  headertext = ...
  ['#!/bin/sh \n' ...
  '# Your job name \n' ...
  '#$ -N ',expname,' \n' ...
  '# \n' ...
  '# Use current working directory \n' ...
  '#$ -cwd \n' ...
  '# \n' ...
  '# Join stdout and stderr \n' ...
  '#$ -j y \n' ...
  '# \n' ...
  '# \n' ...
  '# pe request for MPI. Set your number of processors here. \n' ...
  '#$ -pe dc* ',num2str(nodes),' \n' ...
  '# \n' ...
  '# Run job through bash shell \n' ...
  '#$ -S /bin/bash \n' ...
  '# \n' ...
  '## Output file \n' ...
  '#$ -o ./output.txt \n' ...
  '# \n' ...
  '# The following is for reporting only. It is not really needed \n' ...
  '# to run the job. It will show up in your output file. \n' ...
  '# echo "Got $num_proc processors." \n' ...
  '# echo "Machines:" \n' ...
  '# cat $TMPDIR/machines \n' ...
  '# \n' ...
  '# Use full pathname to make sure we are using the right mpirun \n' ...
  '# Gridengine will set NSLOTS and TMPDIR for you \n' ...        
  '#$ -m bea \n' ...
  '#$ -l h_data=3G,h_rt=336:00:00,highp,exclusive \n' ...
  '. /u/local/Modules/default/init/modules.sh \n' ...
  'module load intel/2020.4 \n' ...
  'mpirun ./mitgcmuv'];
  
  %%% Open output script and write header text
  wfid = fopen(fullfile(dirname,'run_mitgcm'),'w');
  if (wfid == -1)
    error('Could not open PBS script file');
  end
  fprintf(wfid,headertext); 
  
  %%% Close files when we're done

end