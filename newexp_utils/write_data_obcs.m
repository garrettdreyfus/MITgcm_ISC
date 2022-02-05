%%%
%%% write_data_obcs
%%%
%%% Writes the 'data.obcs' input file from 
%%% the cell array OBCS_PARM of parmlist objects
%%% 'realfmt' specifies the output format for real values,
%%% and 'listterm' specifies the fortran NAMELIST termination character.
%%%
function write_data_obcs (dirname,OBCS_PARM,listterm,realfmt)

  %%% Open the 'data.obcs' file for writing
  fname = 'data.obcs';
  fid = fopen(fullfile(dirname,fname),'w');
  if (fid == -1)
    error(['Could not open ',fname,' file']);
  end
  
  %%% Write out file title text
  titletext = [
    '# ===================\r\n' ...
    '# | OBCS parameters |\r\n' ...
    '# ===================\r\n' ...
    '\r\n' ];
  fprintf(fid,titletext);
  
    %%% Parameter section titles
  titles = {...
    'Open boundary types and locations', ...
    'Orlanski parameters', ...
    'Sponge-layer parameters', ...
    'Stevens parameters', ...
    'Sea ice sponge parameters'};
  
    %%% For each parameter section 
  for nparam=1:1:length(OBCS_PARM)
    
    %%% Write section header
    fprintf(fid,['# ',titles{nparam},'\r\n']);
    fprintf(fid,[' &OBCS_PARM0',num2str(nparam),'\r\n']);   
        
    %%% Write each parameter out to the 'data' file
    nextparmlist = OBCS_PARM{nparam};
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


