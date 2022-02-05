%%%
%%% writeDatasetA.m
%%%
%%% Writes an input dataset for MITgcm to a binary input file. This version
%%% of the function appends data to a previous file.
%%%
function writeDatasetA (data,filename,format,precision,append)

  if (append)
    fid=fopen(filename,'a',format); 
  else
    fid=fopen(filename,'w',format); 
  end
  if (fid == -1)
    error(['Could not open ',filename]);
  end
  fwrite(fid,data,precision); 
  fclose(fid);

end

