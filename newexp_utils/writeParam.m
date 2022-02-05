%%%
%%% writeParam.m
%%%
%%% Writes the parameter object parmobj to the file specified  by
%%% 'fid'. 'realfmt' specifies the formatting for real numbers.
%%%
function writeParam (fid,parmobj,realfmt)

  %%% Get parameter type definitions
  paramTypes;
  
  %%% Load parameter properties
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
      paramStr = ['''',paramVal,''''];    
    case PARM_BOOL
      paramStr = inputBool(paramVal);
    case PARM_INTS
      paramStr = list2str(paramVal,'%d');
    case PARM_REALS
      paramStr = list2str(paramVal,realfmt);
    otherwise
      error(['Parameter type ',paramType,' not recognised for parameter ',paramName]);
  end      
  fprintf(fid,[' ',paramName,'=',paramStr,',\r\n']);
      
end

