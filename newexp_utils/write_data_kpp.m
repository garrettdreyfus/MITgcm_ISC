%%%
%%% write_data_kpp
%%%
%%% Writes the 'data.kpp' input file from 
%%% the cell array KPP_PARM of parmlist objects
%%% 'realfmt' specifies the output format for real values,
%%% and 'listterm' specifies the fortran NAMELIST termination character.
%%%
function write_data_kpp (dirname,KPP_PARM,listterm,realfmt)

  %%% Open the 'data.kpp' file for writing
  fname = 'data.kpp';
  fid = fopen(fullfile(dirname,fname),'w');
  if (fid == -1)
    error(['Could not open ',fname,' file']);
  end
  
  %%% Write out file title text
  titletext = [
    '# ===================\r\n' ...
    '# | KPP parameters  |\r\n' ...
    '# ===================\r\n' ...
    '\r\n' ];
  fprintf(fid,titletext);
  
    %%% Parameter section titles
  titles = {...
    'I/O related parameters', ...
    'Genral KPP parameters', ...
    'Boundary layer parameters (S/R bldepth)', ...
    'Boundary layer mixing parameters (S/R blmix)',...
    'Interior mixing parameters (S/R Ri_iwmix)',...
    'Double-diffusive mixing parameters (S/R ddmix)'};
  
    %%% For each parameter section 
  for nparam=1:1:length(KPP_PARM)
    
    %%% Write section header
    fprintf(fid,['# ',titles{nparam},'\r\n']);
    fprintf(fid,[' &KPP_PARM0',num2str(nparam),'\r\n']);   
        
    %%% Write each parameter out to the 'data' file
    nextparmlist = KPP_PARM{nparam};
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