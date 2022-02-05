%%%
%%% inputBool.m
%%%
%%% Formats a boolean as an MITgcm input parameter.
%%%
function strbool = inputBool( bool )

  if (bool)
    strbool = '.TRUE.';
  else
    strbool = '.FALSE.';
  end

end

