%%%
%%% createPBSfile.m
%%%
%%% Writes out a script that runs the MITgcm. Takes as arguments the name
%%% of the experiment, the number of nodes that will be used, and the wall
%%% time in hours (after which the computation will be cut off regardless 
%%% of whether it is complete).
%%%
function createPBSfile_Gordon (dirname,expname,nodes,walltime,account,cluster_path)

  %%% Limit walltime
  walltime = min(ceil(walltime),48);
  
  %%% Set variables correctly for PBS
  headertext = [...
  '#!/bin/bash \n' ...
  '#PBS -q normal \n' ...
  '#PBS -l nodes=',num2str(ceil(nodes/16)),':ppn=16:native \n' ...
  '#PBS -l walltime=',num2str(walltime),':00:00 \n' ...
  '#PBS -N ',expname,' \n' ...
  '#PBS -o output.txt \n' ...
  '#PBS -e errors.txt \n' ...
  '#PBS -A ',account,' \n' ...
  '#PBS -M jhazel@atmos.ucla.edu \n' ...
  '#PBS -m abe \n' ...
  '#PBS -V \n' ...
  '\n' ...
  'cd ',cluster_path,'\n' ...
  'mpirun_rsh -np ',num2str(nodes),' -hostfile $PBS_NODEFILE ./mitgcmuv'];
  
  %%% Open template script
  templatename = './DEFAULTS/results/run_mitgcm';
  tfid = fopen(templatename,'r');
  if (tfid == -1)
    error(['Could not read template PBS run file: ',templatename]);
  end
  
  %%% Open output script and write header text
  wfid = fopen(fullfile(dirname,'run_mitgcm'),'w');
  if (wfid == -1)
    error('Could not open PBS script file');
  end
  fprintf(wfid,headertext); 
  
  %%% Copy over all template text
  count = 0;
  while true
    count = count+1;
    nextline = fgetl(tfid);    
    if ~ischar(nextline)
      break;
    end
    if (count <= 5) %%% skip placeholder comments at start of template
      continue;
    end
    fprintf(wfid,'%s\n',nextline);    
  end
  
  %%% Close files when we're done
  fclose(wfid);
  fclose(tfid);

end

