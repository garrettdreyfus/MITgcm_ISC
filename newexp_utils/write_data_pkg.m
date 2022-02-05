%%%
%%% write_data_pkg
%%%
%%% Writes the 'data.pkg' input file from
%%% the cell array LAYERS_PARM of parmlist objects.
%%% 'realfmt' specifies the output format for real values,
%%% and 'listterm' specifies the fortran NAMELIST termination character.
%%%
function write_data_pkg (dirname,PACKAGES_PARM,listterm,realfmt)

  %%% Open the 'data.pkg' file for writing
  fname = 'data.pkg';
  fid = fopen(fullfile(dirname,fname),'w');
  if (fid == -1)
    error(['Could not open ',fname,' file']);
  end
  
  %%% Write out file title text
  titletext = [
    '# =======================\r\n' ...
    '# | Package parameters |\r\n' ...
    '# =======================\r\n' ...
    '\r\n' ];
  fprintf(fid,titletext);
  
  %%% Parameter section titles
  titles = {'Active packages'};
  
  %%% Write section header
  fprintf(fid,['# ',titles{1},'\r\n']);
  fprintf(fid,[' &PACKAGES','\r\n']);   

  %%% Write each parameter out to the 'data' file
  nextparmlist = PACKAGES_PARM{1};
  for n=1:1:nextparmlist.getLength()      
    writeParam(fid,nextparmlist.getParm(n),realfmt);
  end    

  %%% Write section footer
  fprintf(fid,[' ',listterm,'\r\n']);
  fprintf(fid,'\r\n');    
  
  %%% Close the file when we're finished
  fclose(fid);

end


