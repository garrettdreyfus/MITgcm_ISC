%%%
%%% rdmdsWrapper.m
%%%
%%% Convenience wrapper for rdmds.m that tries to read the specified file
%%% 'fname' and iteration 'dumpIter'. If that name/iteration can't be
%%% found, we try dumpIter+1 and dumpIter-1 in case of rounding errors
%%% between Matlab and MITgcm.
%%%
function A = rdmdsWrapper(fname,dumpIter)
  
  A = rdmds(fname,dumpIter);        
  if (isempty(A))
    A = rdmds(fname,dumpIter-1);      
  end
  if (isempty(A))
    A = rdmds(fname,dumpIter+1);      
  end
  
end

