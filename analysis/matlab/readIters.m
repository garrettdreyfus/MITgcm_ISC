%%%
%%% readIters.m
%%%
%%% Reads and sums all iterations of a specified MITgcm output field
%%% between specified times, and calculates the time average.
%%%
function avg = readIters (exppath,field,dumpIters,nDumps,deltaT,tmin,tmax,Nx,Ny,Nr)
 
  avg = zeros(Nx,Ny,Nr);
  navg = 0;
  
  %%% Loop through output iterations
  for n=1:length(dumpIters)
     
    tyears =  dumpIters(n)*deltaT/86400/365;
    
    if ((tyears >= tmin) && (tyears <= tmax))
  
      %%% Read the next iteration and check the data were found      
      avg_tmp = rdmdsWrapper(fullfile(exppath,'results',field),dumpIters(n));            
      if (isempty(avg_tmp))
        disp(['Ran out of data at ITER=',dumpIters(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.'])
        break;
      end
      
      avg = avg + avg_tmp;
      navg = navg + 1;

    end
    
  end
  
  %%% Calculate average
  if (navg > 0)
    avg = avg / navg;
  else
    error('No output files found');
  end

end

