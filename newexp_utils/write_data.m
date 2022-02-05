%%%
%%% write_data
%%%
%%% Writes the 'data' input file from the cell array PARM of parmlist 
%%% objects. 'realfmt' specifies the output formatting for real 
%%% numbers, and 'listterm' specifies the termination character for fortran
%%% NAMELISTS.
%%%
function write_data (dirname,PARM,listterm,realfmt)
  
  %%% Open the 'data' file
  fname = 'data';
  fid = fopen(fullfile(dirname,fname),'w');
  if (fid == -1)
    error(['Could not open ',fname,' file']);
  end
  
  %%% Write out file title text
  titletext = [
    '# ====================\r\n' ...
    '# | Model parameters |\r\n' ...
    '# ====================\r\n' ...
    '\r\n'];
  fprintf(fid,titletext);
  
  %%% Parameter section titles
  titles = {...
    'Continuous equation parameters', ...
    'Elliptic solver parameters', ...
    'Time stepping parameters', ...
    'Gridding parameters', ...
    'Input datasets'};
    
  %%% For each parameter section 
  for nparam=1:1:length(PARM)
    
    %%% Write section header
    fprintf(fid,['# ',titles{nparam},'\r\n']);
    fprintf(fid,[' &PARM0',num2str(nparam),'\r\n']);   
        
    %%% Write each parameter out to the 'data' file
    nextparmlist = PARM{nparam};
    for n=1:1:nextparmlist.getLength()      
      writeParam(fid,nextparmlist.getParm(n),realfmt);
    end
    
    %%% Write section footer
    fprintf(fid,[' ',listterm,'\r\n']);
    fprintf(fid,'\r\n');
    
  end
  
  %%% Close the file when we're finished
  fclose(fid);

end

