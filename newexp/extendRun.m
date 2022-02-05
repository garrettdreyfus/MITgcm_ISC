%%%
%%% extendRun.m
%%%
%%% Extends an MITgcm simulation by restarting from the last checkpoint for a selected period of time.
%%%
%%% expdir - full path to folder containing the experiment
%%% expname - name of the experiment
%%% exttime - new end time for the simulation, in seconds
%%% doRestart - set true to actually execute the simulation in this script
%%%
function extendRun (expdir,expname,newEndTime,doRestart)
 
  %%% Directories containing simulation files  
  inputdir = fullfile(expdir,expname,'input');  
  resultsdir = fullfile(expdir,expname,'results');  
  datafile = fullfile(inputdir,'data');
  paramsfile = fullfile(inputdir,'params.m');
  tmpparamsfile = fullfile(inputdir,'tmp_params.m');
  
  %%% Find the last pickup file in the directory
  resultsFiles = dir(resultsdir);
  newStartIter = -1;
  for m = 1:length(resultsFiles)
    pname = resultsFiles(m).name;
    if (startsWith(pname,'pickup.'))
      startIter = str2num(pname(8:17));
      if (startIter > newStartIter)
        newStartIter = startIter;
      end
    end
  end
 
  %%% Set nIter0 to timestepnumber in the simulation's input data file
  fid = fopen(datafile,'r');
  if (fid == -1)
    error(['Could not open ',datafile]);
  end
  datastr = '';
  tline = fgetl(fid);
  while (ischar(tline))
    if (~isempty(strfind(tline,'nIter0')))
      datastr = [datastr,' nIter0=',num2str(newStartIter),',\n'];
    else
      if (~isempty(strfind(tline,'endTime')))
        datastr = [datastr,' endTime=',num2str(newEndTime),',\n'];
      else
        datastr = [datastr,tline,'\n'];
      end
    end      
    tline = fgetl(fid);
  end
  fclose(fid);
  fid = fopen(datafile,'w');
  fprintf(fid,datastr);
  fclose(fid);
  
  %%% Update the simulation's params.m file to reflect the new simulation
  %%% end time
  fid = fopen(paramsfile,'r');  
  if (fid == -1)
    error(['Could not open ',paramsfile]);
  end
  tfid = fopen(tmpparamsfile,'w');
  if (tfid == -1)
    error(['Could not open ',tmpparamsfile]);
  end
  tline = fgetl(fid);
  while (ischar(tline))
    if (~isempty(strfind(tline,'endTime')))
      paramsstr = ['endTime=',num2str(newEndTime),';'];
    else
      paramsstr = tline;
    end
    fprintf(tfid,'%s\n',paramsstr);
    tline = fgetl(fid);
  end
  fclose(fid);
  fclose(tfid);  
  copyfile(tmpparamsfile,paramsfile);
  delete(tmpparamsfile);
  
  %%% Restart simulation if required
  if (doRestart)
    currentdir = pwd;
    cd(resultsdir);
    system('sh run.sh');
    cd(currentdir);
  end
    
end


