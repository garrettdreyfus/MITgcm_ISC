%%%
%%% parm.m
%%%
%%% Object representing a single MITgcm parameter.
%%%
classdef parm < handle
  
  properties (GetAccess=private)
    name;
    value;
    type;
  end
  
  methods    
    function obj = parm (name, value, type)
      obj.name = name;
      obj.value = value;
      obj.type = type;
    end
    
    function name = getName (obj)
      name = obj.name;
    end
    
    function value = getValue (obj)
      value = obj.value;
    end
    
    function type = getType (obj)
      type = obj.type;
    end    
    
  end
  
end

