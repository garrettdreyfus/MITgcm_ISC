%%%
%%% writeDataset.m
%%%
%%% Writes an input dataset for MITgcm to a binary input file.
%%%
function writeDataset (data,filename,format,precision)

  fid=fopen(filename,'w',format); 
  if (fid == -1)
    error(['Could not open ',filename]);
  end
  fwrite(fid,data,precision); 
  fclose(fid);

end

