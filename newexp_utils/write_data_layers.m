%%%
%%% write_data_layers
%%%
%%% Writes the 'data.layers' input file from
%%% the cell array LAYERS_PARM of parmlist objects.
%%% 'realfmt' specifies the output format for real values,
%%% and 'listterm' specifies the fortran NAMELIST termination character.
%%%
function write_data_layers (dirname,LAYERS_PARM,listterm,realfmt)

  %%% Open the 'data.layers' file for writing
  fname = 'data.layers';
  fid = fopen(fullfile(dirname,fname),'w');
  if (fid == -1)
    error(['Could not open ',fname,' file']);
  end
  
  %%% Write out file title text
  titletext = [
    '# =====================\r\n' ...
    '# | Layers parameters |\r\n' ...
    '# =====================\r\n' ...
    '\r\n' ];
  fprintf(fid,titletext);
  
  %%% Parameter section titles
  titles = {'Parameters for layers package'};
  
  %%% For each parameter section 
  for nparam=1:1:length(LAYERS_PARM)
    
    %%% Write section header
    fprintf(fid,['# ',titles{nparam},'\r\n']);
    fprintf(fid,[' &LAYERS_PARM0',num2str(nparam),'\r\n']);   
        
    %%% Write each parameter out to the 'data' file
    nextparmlist = LAYERS_PARM{nparam};
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


