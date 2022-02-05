%%%
%%% extendBatch.m
%%% 
%%% Extends a batch of simulations.
%%%
 
%%% Directory containing experiments
% expdir = '/Volumes/LaCie/UCLA/Projects/MITgcm_ACC_AABW/experiments';
expdir = '/scratch/03198/astewart/MITgcm_ACC_AABW/experiments';
 
%%% List of experiments to restart
expnames = { ...  
  'ACC_AABW_lores_Hr250_Ht0_kap1e-3', ...
  'ACC_AABW_lores_Hr750_Ht0_kap1e-3', ... 
  'ACC_AABW_lores_Hr1250_Ht0_kap1e-3', ... 
  'ACC_AABW_lores_Hr1000_Ht0_kap1e-3_taue0.1', ... 
  'ACC_AABW_lores_Hr1000_Ht0_kap1e-3_taue0.1', ... 
  'ACC_AABW_lores_Hr1000_Ht0_kap1e-3_Taabw-0.75', ... 
  'ACC_AABW_lores_Hr1000_Ht0_kap1e-3_Taabw-0.5', ... 
};
  
%%% Temporal parameters
t1year = 86400*365;
newEndTime = 40*t1year;
 
%%% Set true to actually restart the simulations now
doRestart = false;
 
%%% Restart all simulations
for n = 1:length(expnames)
  extendRun(expdir,expnames{n},newEndTime,doRestart);
end


