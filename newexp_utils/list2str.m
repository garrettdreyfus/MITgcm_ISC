%%%
%%% list2str.m
%%%
%%% Formats a list of numbers as an input parameter for MITgcm.
%%%
function str = list2str (list,fmt)

  %%% Not the most efficient way of doing this, but generally this
  %%% operation isn't performance-limiting 
  lfreq = 10;
  str = '';
  for n=1:1:length(list)      
    str = [str,num2str(list(n),fmt),];
    if (n ~= length(list))
      str = [str,', '];
    end
    if ((mod(n,lfreq)==0) && (n ~=length(list)))
      str = [str,'\r\n      '];
    end
  end

end

