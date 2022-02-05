%%%
%%% write_data_diagnostics.m
%%%
%%% Writes the 'data.diagnostics' input file from
%%% the cell array DIAG_PARM of parameter list objects.
%%% 'realfmt' specifies the output format for real values,
%%% and 'listterm' specifies the fortran NAMELIST termination character.
%%%
function write_data_diagnostics (dirname,DIAG_PARM,listterm,realfmt)

  %%% Open the 'data.diagnostics' file for writing
  fname = 'data.diagnostics';
  fid = fopen(fullfile(dirname,fname),'w');
  if (fid == -1)
    error(['Could not open ',fname,' file']);
  end
  
  %%% Write out file title text
  titletext = [
    '# ==========================\r\n' ...
    '# | DIAGNOSTICS parameters |\r\n' ...
    '# ==========================\r\n' ...
    '\r\n' ];
  fprintf(fid,titletext);

  for nparam=1:2

    %%% Write section header
    if (nparam==1)
      fprintf(fid,['# Diagnostics package choices\n']);
      fprintf(fid,[' &DIAGNOSTICS_LIST','\n']);   
    else
      fprintf(fid,['# Statistics package choices\n']);
      fprintf(fid,[' &DIAG_STATIS_PARMS','\n']);
    end

    %%% Write each parameter out to the 'data' file    
    nextparmlist = DIAG_PARM{nparam};
    for n=1:1:nextparmlist.getLength()      
      writeParam(fid,nextparmlist.getParm(n),realfmt);
    end    

    %%% Write section footer
    fprintf(fid,[' ',listterm,'\n']);
    fprintf(fid,'\r\n');
    
  end
  
  %%% Close the file when we're finished
  fclose(fid);

end


