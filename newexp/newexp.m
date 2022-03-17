%%%
%%% newexp.m
%%%
%%% Convenience script to set up folders and options for a new MITgcm run
%%% in a way that allows:
%%% - automatic generation of basic folder and file structure
%%% - automatic calculation of inter-dependent input parameters
%%% - compatibility between required numbers of processors
%%%
%%% Sets up a new experiment folder with subdirectories for the build, the 
%%% code opti ons, the input files, and a results folder. Creates a 'build'
%%% script in the build folder and a 'run' script in the results folder.
%%% Generates the SIZE.h file in the code folder based on parameters
%%% specified here, and copies other code files from the DEFAULTS/code 
%%% folder. Generates all 'eedata' and some 'data' parameters in the input 
%%% folder based on parameters specified here and code in create_data 
%%% function. Other parameters are copied from the DEFAULTS/input folder.
%%%
%%% NOTE: 'expname' MUST NOT be set to 'DEFAULTS'
%%%

 function newexp (batch_name,exp_name)

   
%    addpath /data/MITgcm_ASF-csi/newexp_utils/;
  addpath ../newexp_utils/;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% USER-SPECIFIED OPTIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%clo%%%%%%%%%%  
  
  
  %%% Local directory in which to create experiments
  runsdir = '/Volumes/Kilchoman/UCLA/Projects/MITgcm_ISC/';
 
  %%% Experiment subdirectories 
  builddir = 'build';
  codedir = 'code';
  inputdir = 'input';
  resultsdir = 'results';  
 
  %%% List terminator character for parameter files - may be '/' or '&'
  %%% depending on operating system
  listterm = '&';

  %%% Line feed character - important for .sh shell script files
  %%% On unix this probably need to be '\n', in windows '\r\n\'
  lf = '\n';    
  
  nSx = 1; %%% no. of tiles per processor in x-direction
  nSy = 1; %%% no. of tiles per processor in y-direction
  nTx = 1; %%% no. of threads per processor in x-direction
  nTy = 1; %%% no. of threads per processor in y-direction
  OLx = 3; %%% no. of overlapping x-gridpoints per tile
  OLy = 3; %%% no. of overlapping y-gridpoints per tile    
  

  
  %%% These parameters are most likely to vary between experiments
  %%% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  
  
  %%% Set-up for Ardbeg  
%   opt_file = 'ardbeg_ucla'; %%% options file name
%   use_mpi = true; %%% set true for parallel processing
%   use_pbs = true; %%% set true for execution via PBS
%   cluster = 'ardbeg';
%   queue = 'all.q'; 
%   sNx = 20; %%% no. of x-gridpoints per tile
%   sNy = 25; %%% no. of y-gridpoints per tile
%   nPx = 10; %%% no. of processors in x-direction
%   nPy = 8; %%% no. of processors in y-direction
%   Nr = 100; %%% no. of z-gridpoints

  %%% Set-up for Hoffman2
  opt_file = 'hoffman2_ucla'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  cluster = 'hoffman2';
  queue = 'all.q';
  sNx = 20; %%% no. of x-gridpoints per tile
  sNy = 25; %%% no. of y-gridpoints per tile
  nPx = 10; %%% no. of processors in x-direction
  nPy = 8; %%% no. of processors in y-direction
  Nr = 100; %%% no. of z-gridpoints
  
  
  %%% N.B. bverything in this section below this is old code, and will not
  %%% work
  
%   %%% Set-up for Si's Mac Pro - barotropic test case
%   opt_file = 'darwin_amd64_gfortran'; %%% options file name
%   use_mpi = false; %%% set true for parallel processing
%   use_pbs = false; %%% set true for execution via PBS
%   cluster = 'none';    
%   sNx = 300; %%% no. of x-gridpoints per tile
%   sNy = 250; %%% no. of y-gridpoints per tile
%   nPx = 1; %%% no. of processors in x-direction
%   nPy = 1; %%% no. of processors in y-direction
%   Nr = 40; %%% no. of z-gridpoints
    
  
% %%% Set-up for Stampede2
%   opt_file = 'xsede_stampede'; %%% options file name
%   use_mpi = true; %%% set true for parallel processing
%   use_pbs = true; %%% set true for execution via PBS
%   cluster = 'stampede';   
%   queue = 'normal';  
% %   queue = 'flat-quadrant';
% 
%  if (is_hires) 
%   %%% 1 node,400*448
%   sNx = 50; %%% no. of x-gridpoints per tile
%   sNy = 56; %%% no. of y-gridpoints per tile
%   nPx = 8; %%% no. of processors in x-direction
%   nPy = 8; %%% no. of processors in y-direction
%   Nr = 70; %%% no. of z-gridpoint 
% 
% %   sNx = 50; %%% no. of x-gridpoints per tile
% %   sNy = 28; %%% no. of y-gridpoints per tile
% %   nPx = 8; %%% no. of processors in x-direction
% %   nPy = 16; %%% no. of processors in y-direction
% %   Nr = 70; %%% no. of z-gridpoint 
%   
% %   sNx = 25; %%% no. of x-gridpoints per tile
% %   sNy = 28; %%% no. of y-gridpoints per tile
% %   nPx = 16; %%% no. of processors in x-direction
% %   nPy = 16; %%% no. of processors in y-direction
% %   Nr = 70; %%% no. of z-gridpoint 
% 
% %   sNx = 36; %%% no. of x-gridpoints per tile
% %   sNy = 37; %%% no. of y-gridpoints per tile
% %   nPx = 11; %%% no. of processors in x-direction
% %   nPy = 12; %%% no. of processors in y-direction
% 
%  else   
% %%%% ~ 2km grid spacing
%   sNx = 25; %%% no. of x-gridpoints per tile
%   sNy = 28; %%% no. of y-gridpoints per tile
%   nPx = 8; %%% no. of processors in x-direction
%   nPy = 8; %%% no. of processors in y-direction
%   Nr = 70; %%% no. of z-gridpoints
% 
%  end
%  
%   acct = 'tg854737';

 




  %%% Set-up for Gordon
%    opt_file = 'xsede_gordon'; %%% options file name
%    use_mpi = true; %%% set true for parallel processing
%    use_pbs = true; %%% set true for execution via PBS
%    cluster = 'gordon';
%    queue = 'normal';
%    acct = 'ddp258';

 %%% Set-Up for Comet
%    opt_file = 'xsede_comet'; %%% options file name
%    use_mpi = true; %%% set true for parallel processing
%    use_pbs = true; %%% set true for execution via PBS
%    cluster = 'comet';
%    queue = 'normal';
%    acct = 'ddp258';
   

  %%% Uploading/downloading parameters 
  switch(cluster)
            
    case 'gordon'
        
      username = 'csi';
      clustername = 'gordon.sdsc.xsede.org';
      toolsdir = fullfile('/oasis/projects/nsf/',acct,username,'/MITgcm_WS/tools/');
      clusterdir = fullfile('/oasis/projects/nsf/',acct,username,'/MITgcm_WS/',batch_name);      
           
    case 'stampede'
        
      username = 'tg854737';
      clustername = 'stampede2.tacc.utexas.edu';
      toolsdir = fullfile('/scratch/06174/tg854737/MITgcm_ASF-csi/tools/');
      clusterdir = fullfile('/scratch/06174/tg854737/MITgcm_ASF-csi/',batch_name);      
      
    case 'comet'
        
      username = 'csi';
      clustername = 'comet.sdsc.xsede.org';
      toolsdir = fullfile('/oasis/projects/nsf/',acct,username,'/MITgcm_WS/tools/');
      clusterdir = fullfile('/oasis/projects/nsf/',acct,username,'/MITgcm_WS/experiments/',batch_name);      
            
    case 'ardbeg' 
  
      username = 'astewart';
      clustername = 'caolila.atmos.ucla.edu';
      toolsdir = '/data2/garrett/MITgcm_ISC/tools/';
      clusterdir = fullfile('/data2/garrett/MITgcm_ISC/experiments/',batch_name);      
      
    case 'hoffman2'
        
      username = 'andrewst';
      clustername = 'hoffman2.idre.ucla.edu';
      toolsdir = '/u/scratch/a/andrewst/MITgcm_ISC/tools/';
      clusterdir = ['/u/scratch/a/andrewst/MITgcm_ISC/experiments/',batch_name];
      
    otherwise %%% Defaults to 'none'
       
      username = 'csi';
      clustername = '';
      toolsdir = '/Users/csi/MITgcm_ASF-csi/tools/';
      clusterdir = fullfile('/Users/csi/MITgcm_ASF-csi/experiments/',batch_name);  
     
      
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% NON-USER-SPECIFIED PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Paths to sub-directories
  dirpath = fullfile(runsdir,batch_name);
  exppath = fullfile(dirpath,exp_name);
  buildpath = fullfile(exppath,builddir);
  codepath =  fullfile(exppath,codedir);
  inputpath = fullfile(exppath,inputdir);
  resultspath = fullfile(exppath,resultsdir);
  
  %%% We have to use MPI if we're using PBS
  if (use_pbs)
    use_mpi = true;    
  end
  
  %%% If use_mpi is false then we can only have one processor
  if ((use_mpi == false) && ((nPx~=1) || (nPy~=1)))
    error('Only one processor allowed for non-MPI computation');
  end
  
  %%% Calculate total grid size and number of nodes
  Nx = sNx*nSx*nPx;
  Ny = sNy*nSy*nPy;
  nodes = nPx*nPy;



  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIRECTORIES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%


  %%% Open experiment top-level directory
  [dir_success,dir_msg,dir_msgid] = mkdir(exppath);
  if (dir_success == 0)
    error(strcat(['Could not open ',exp_name,' : ',num2str(dir_msgid),' : ',dir_msg]));
  end

  %%% Open sub-directories
  subdirnames = {builddir,codedir,inputdir,resultsdir};
  for n=1:1:length(subdirnames)     
    [dir_success,dir_msg,dir_msgid] = mkdir(exppath,subdirnames{n});
    if (dir_success == 0)
      error(strcat(['Could not open ',exppath,subdirnames{n},' : ',num2str(dir_msgid),' : ',dir_msg]));
    end
  end


  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INPUT FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%

  
  %%% Copy other files across
  codelist = dir('./DEFAULTS/input/');
  for n=1:1:length(codelist)
    %%% Ignore hidden files
    if (codelist(n).name(1) == '.')
      continue;
    end    
    copyfile(fullfile('./DEFAULTS/input/',codelist(n).name),fullfile(inputpath,codelist(n).name));
  end   

  %%% Generate experiment input files ('data','data.seaice', etc.)
  nTimeSteps = setParams(exp_name,inputpath,codepath,listterm,Nx,Ny,Nr);  

  %%% Generate 'eedata'
  create_eedata(inputpath,listterm,nTx,nTy);     

  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% CODE FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%


  %%% Generate SIZE.h and just copy other code files
  createSIZEh(codepath,sNx,sNy,nSx,nSy,nPx,nPy,OLx,OLy,Nr);
  codelist = dir('./DEFAULTS/code/');
  for n=1:1:length(codelist)
    if (codelist(n).name(1) == '.')
      continue;
    end
    copyfile(fullfile('./DEFAULTS/code/',codelist(n).name),fullfile(codepath,codelist(n).name));
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ESTIMATE WALL TIME %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Computation time (in hours) per gridpoint (in space and time) 
  %%% assigned to each processor.  
  switch (cluster)
    case 'gordon'
      alpha = 0.63e-9;
    case 'stampede'
      alpha = 1.59e-9;
    case 'comet'
      alpha = 0.63e-9;
    case 'hoffman2'
      alpha = 0.63e-9;
    otherwise %%% Defaults to Ardbeg
      alpha = 0.5e-9;
  end  
  
  %%% Estimated total computation time in hours (assumes one tile per
  %%% processor). 
  
  %%% This tends to overestimate the computation time when OLx is
  %%% comparable to sNx or OLy is comparable to sNy.
  % comptime = alpha*(sNx+2*OLx)*(sNy+2*OLy)*Nr*Nt  
  
  %%% This seems to provide a decent estimate when OLx is
  %%% comparable to sNx or OLy is comparable to sNy; 'ghost' gridpoints
  %%% require perhaps half as much processing as 'real' gridpoints.
  comptime = alpha*(sNx+OLx)*(sNy+OLy)*Nr*nTimeSteps
    
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR BUILDING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Build commands - depend on whether MPI is used
  if (use_mpi)
    mpistr = '-mpi ';
  else 
    mpistr = '';
  end
  buildcommands = strcat([...    
    'rm my_opt_file',lf,...
    'ln -s ',fullfile(toolsdir,'build_options',opt_file),' my_opt_file',lf, ...
    fullfile(toolsdir,'genmake2'),' ',mpistr,'-mods=../code -of=my_opt_file',lf, ...
    'make depend',lf, ...
    'make --always-make -j 2',lf,]);

  %%% Create the 'build' shell file
  fid = fopen(fullfile(buildpath,'build.sh'),'w');
  if (fid == -1)
    error('Could not open build script for writing');
  end
  fprintf(fid,buildcommands);
  fclose(fid);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR RUNNING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cl%%%%%%%%%%%


  %%% Commands to link input and build folders to results folder
  runcommands = [...
    'ln -s ../',inputdir,'/* ./ ',lf, ...
    'ln -s ../',builddir,'/mitgcmuv ',lf];

  %%% Execution command depends on whether MPI is used
  if (use_pbs)
    switch (cluster)
      case 'gordon'
        createPBSfile_Gordon(resultspath,exp_name,nodes,2*comptime,acct,fullfile(clusterdir,exp_name,resultsdir));        
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
      case 'stampede'
        createPBSfile_Stampede(resultspath,exp_name,nodes,2*comptime,acct,fullfile(clusterdir,exp_name,resultsdir));       
        runcommands = [runcommands,'sbatch run_mitgcm',lf];
      case 'stampede2'
        createPBSfile_Stampede2(resultspath,exp_name,nodes,2*comptime,acct,fullfile(clusterdir,exp_name,resultsdir));       
        runcommands = [runcommands,'sbatch run_mitgcm',lf];
      case 'comet'
        createPBSfile_Comet(resultspath,exp_name,nodes,2*comptime,acct,fullfile(clusterdir,exp_name,resultsdir));        
        runcommands = [runcommands,'sbatch run_mitgcm_comet',lf];
      case 'hoffman2'
        createPBSfile_Hoffman(resultspath,exp_name,nodes);        
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
      otherwise %%% Defaults to Ardbeg
        createPBSfile(resultspath,exp_name,nodes,queue);
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
    end    
  else
    if (use_mpi)
      runcommands = [runcommands,'mpirun -np ',num2str(nodes), ...
                      ' ./mitgcmuv > output.txt',lf];
    else
      runcommands = [runcommands,'./mitgcmuv > output.txt',lf];
    end
  end

  %%% Create the 'run' shell file
  fid = fopen(fullfile(resultspath,'run.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,runcommands);
  fclose(fid);
  
  %%% Copy other files across
  resultslist = dir('./DEFAULTS/results/');
  for n=1:1:length(resultslist)
    %%% Ignore hidden files and run script template
    if ((resultslist(n).name(1) == '.') || strcmp(resultslist(n).name,'run_mitgcm'))
      continue;
    end    
    copyfile(fullfile('./DEFAULTS/results/',resultslist(n).name),fullfile(resultspath,resultslist(n).name));
  end    


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR COMPILING AND RUNNING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% A file to build and run MITgcm, for user convenience if they're
  %%% confident that everything is set up correctly

  %%% Commands to link input and build folders to results folder
  commands = [...
    'cd ./',builddir,'/',lf, ...
    'sh build.sh',lf, ...
    'cd ../',resultsdir,'/ ',lf, ...
    'sh run.sh',lf ];

  %%% Create the 'run' shell file
  fid = fopen(fullfile(exppath,'build_and_run.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,commands);
  fclose(fid);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILES FOR UPLOADING AND DOWNLOADING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Upload command
  upcommand = [...
    'rsync -av --update ',...
    '../',exp_name,' ', ...
    username,'@',clustername,':',clusterdir];
  fid = fopen(fullfile(exppath,'upload_to_cluster.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,upcommand);
  fclose(fid);
  
  %%% Download command
  downcommand = [...
    'rsync -av --update ', ...
    username,'@',clustername,':', ...
    fullfile(clusterdir,exp_name,resultsdir),'/* ', ...
    './results/ \n'];  
  fid = fopen(fullfile(exppath,'download_from_cluster.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,downcommand);  
  fclose(fid);
  
end


