%%%
%%% parmlist.m
%%%
%%% Object representing a grouping of MITgcm parameters, corresponding to 
%%% a Fortran NAMELIST in the MITgcm input.
%%%
classdef parmlist < handle
  
  properties (GetAccess=private)
    parmcell;
    nparms;
  end
  
  methods    
    
    function obj = parmlist ()
      obj.parmcell = {};
      obj.nparms = 0;
    end
    
    function parm = getParm (obj,idx)
      if (idx > obj.nparms)
        parm = [];
      else
        parm = obj.parmcell{idx};
      end
    end
    
    function addParm (obj,name,value,type)
      newparm = parm(name,value,type);
      obj.parmcell = [obj.parmcell,{newparm}];
      obj.nparms = obj.nparms + 1;
    end
    
    function length = getLength (obj)
      length = obj.nparms;
    end
    
  end
  
end

