%%%
%%% avg_t.m
%%%
%%% Calculates the time average of the output fields from MITgcm runs.
%%%

%%% Read experiment data
clear diag_fields;
clear diag_timePhase;
clear diag_fileNames;
clear diag_frequency;
loadexp;

Nlayers = 70;

savename = [prodir '/' expname '_tavg_5yrs.mat'];


%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.-p[      
dumpFreq = diag_frequency(1);
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Calculate time average for whichever fields are present
 
for m=1:length(diag_fields)
% % for m=1:56

    m
  if (m == 1) %%% N.B. This won't work if the first field isn't one of those listed below
    flag = '';
  else    
    flag = '-append';
  end
  var_name = diag_fields{m};
  if  (strcmp(var_name,'LaUH1RHO')|strcmp(var_name,'LaVH1RHO')|...
          strcmp(var_name,'LaHw1RHO')|strcmp(var_name,'LaHs1RHO')|...
          strcmp(var_name,'LaUH1RHO_inst')|strcmp(var_name,'LaVH1RHO_inst')|...
          strcmp(var_name,'LaHw1RHO_inst')|strcmp(var_name,'LaHs1RHO_inst')|...
          strcmp(var_name,'LaUH2TH')|strcmp(var_name,'LaVH2TH')|...
          strcmp(var_name,'LaHw2TH')|strcmp(var_name,'LaHs2TH')|...
          strcmp(var_name,'LaUH2TH_inst')|strcmp(var_name,'LaVH2TH_inst')|...
          strcmp(var_name,'LaHw2TH_inst')|strcmp(var_name,'LaHs2TH_inst'))   
      if (diag_frequency(m) > 0)
          clear avg var_data
        var_data = readIters(exppath,var_name,dumpIters,nDumps,deltaT,tmin,tmax,Nx,Ny,Nlayers);    
        tempStruct.(var_name) = var_data;   
        save(savename,'-struct','tempStruct',var_name,flag);
      end     
  else
      if (diag_frequency(m) > 0)
        var_data = readIters(exppath,var_name,dumpIters,nDumps,deltaT,tmin,tmax,Nx,Ny,Nr);    
        tempStruct.(var_name) = var_data;   
        save(savename,'-struct','tempStruct',var_name,flag);
      end
  end
end




% dumpFreq = diag_frequency(57);
% nDumps = floor(nTimeSteps*deltaT/dumpFreq);
% dumpIters = round((1:nDumps)*dumpFreq/deltaT);
% dumpIters = dumpIters(dumpIters > nIter0);
% 
% 
% for m=57:length(diag_fields)
%     m
%   if (m == 1) %%% N.B. This won't work if the first field isn't one of those listed below
%     flag = '';
%   else    
%     flag = '-append';
%   end
%   var_name = diag_fields{m};
%   if  (strcmp(var_name,'LaUH1RHO')|strcmp(var_name,'LaVH1RHO')|...
%           strcmp(var_name,'LaHw1RHO')|strcmp(var_name,'LaHs1RHO')|...
%           strcmp(var_name,'LaUH1RHO_inst')|strcmp(var_name,'LaVH1RHO_inst')|...
%           strcmp(var_name,'LaHw1RHO_inst')|strcmp(var_name,'LaHs1RHO_inst')|...
%           strcmp(var_name,'LaUH2TH')|strcmp(var_name,'LaVH2TH')|...
%           strcmp(var_name,'LaHw2TH')|strcmp(var_name,'LaHs2TH')|...
%           strcmp(var_name,'LaUH2TH_inst')|strcmp(var_name,'LaVH2TH_inst')|...
%           strcmp(var_name,'LaHw2TH_inst')|strcmp(var_name,'LaHs2TH_inst'))   
%       if (diag_frequency(m) > 0)
%           clear avg var_data
%         var_data = readIters(exppath,var_name,dumpIters,deltaT,tmin,tmax,Nx,Ny,Nlayers);    
%         tempStruct.(var_name) = var_data;   
%         save(savename,'-struct','tempStruct',var_name,flag);
%       end     
%   else
%       if (diag_frequency(m) > 0)
%         var_data = readIters(exppath,var_name,dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);    
%         tempStruct.(var_name) = var_data;   
%         save(savename,'-struct','tempStruct',var_name,flag);
%       end
%   end
% end
