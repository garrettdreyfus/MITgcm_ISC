%%%
%%% write_matlab_params.m
%%%
%%% Writes out the input parameter list for the MITgcm run as Matlab script
%%% file, for ease of post-processing. PARM must be a cell array of
%%% parmlist objects.
%%%
function write_matlab_params(dirname,PARM,realfmt)
  
  %%% Get parameter type definitions
  paramTypes;
  
  %%% Open the file
  fname = 'params.m';
  fpath = fullfile(dirname,fname);
  fid = fopen(fpath,'w');
  if (fid == -1)
    error(['Could not open ',fpath]);
  end
  
  %%% Write documentation
  headertext = [...
    '%%%%%% \r\n', ...
    '%%%%%% ',fname,'\r\n' ...
    '%%%%%% \r\n' ...
    '%%%%%% Automatically-generated Matlab file containing parameters\r\n' ...
    '%%%%%% used in this MITgcm experiment.\r\n' ...
    '%%%%%% \r\n' ...
    '\r\n'];
  fprintf(fid,headertext);
  
  %%% For each parameter section 
  for nparam=1:1:length(PARM)                
    
    parmlistobj = PARM{nparam};
    listlen = getLength(parmlistobj);
    
    %%% Write each parameter out to the 'data' file
    for n=1:listlen   
  
      parmobj = parmlistobj.getParm(n);
      paramName = parmobj.getName();
      paramVal = parmobj.getValue();
      paramType = parmobj.getType();
      
      %%% Write the parameter based on its type
      switch paramType;
        case PARM_INT
          paramStr = num2str(paramVal);
        case PARM_REAL
          paramStr = num2str(paramVal,realfmt);
        case PARM_STR          
          paramName(paramName=='(') = '{';
          paramName(paramName==')') = '}';          
          paramStr = ['''',paramVal,''''];        
        case PARM_BOOL
          paramStr = num2str(paramVal~=0);
        case PARM_INTS
          paramStr = ['[ ',num2str(paramVal,'%d '),' ]'];
        case PARM_REALS          
          paramStr = ['[ ',num2str(reshape(paramVal,[1 length(paramVal)]),[realfmt,' ']),' ]'];
        otherwise
          error(['Parameter type ',paramType,' not recognised for parameter ',paramName]);
      end    
      fprintf(fid,[paramName,'=',paramStr,';\r\n']);
      
    end
    
    %%% Write section footer    
    fprintf(fid,'\r\n');
    
  end
  
  %%% Close the output file
  fclose(fid);

end

