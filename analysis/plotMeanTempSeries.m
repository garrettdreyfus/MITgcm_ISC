%%%
%%% plotMeanTempSeries.m
%%%
%%% Plots the instantaneous domain-mean potential temperature from 
%%% MITgcm simulations.
%%%
%%%

%%% Read experiment data
loadexp;

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(end);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters >= nIter0);
nDumps = length(dumpIters);

%%% Matrices for integration
DX = repmat(delX',[1 Ny Nr]);
DY = repmat(delY,[Nx 1 Nr]);
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

tt = zeros(1,nDumps);
theta_avg = zeros(1,nDumps);
ptlen = 0;
for n=1:nDumps
% for n=364:364
 
  tt(n) =  dumpIters(n)*deltaT/86400/365;
  
  %%% Attempt to load either instantaneous velocities or their squares
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA_inst'),dumpIters(n)) ;      
 
  if (isempty(theta))   
    break;
  end  
  
  %%% Calculate domain-mean potential temperature
  theta_avg(n) = sum(sum(sum(theta.*DX.*DY.*DZ.*hFacC))) / sum(sum(sum(DX.*DY.*DZ.*hFacC)));
  
  %%% Increment counter
  ptlen = ptlen + 1;
  
end

disp(expname);
disp(tt(n));

%%% Calculate trend
trend_length = 20;
if (ptlen >= trend_length)
  p = polyfit(tt(ptlen-trend_length+1:ptlen),theta_avg(ptlen-trend_length+1:ptlen),1);
  disp(['Trend in ',expname,' = ',num2str(p(1)*100),' deg C/century']);
end

if (~isempty(find(isnan(theta_avg))))
  disp(['NaNs found in experiment ',expname]);
end

figure(1);
clf;
axes('FontSize',16);
plot(tt(1:ptlen),theta_avg(1:ptlen),'o-');
plot(tt(1:ptlen),theta_avg(1:ptlen)-theta_avg(1),'o-');
axis tight;
xlabel('t (years)');
ylabel('Mean potential temperature (^oC)');

