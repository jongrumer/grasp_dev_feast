%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
     
function M = nuclear(mu,Q,I)
     
% This function computes the reduced nuclear matrix elements
% in atomic units. 
%
% Written by Per Jonsson and Martin Andersson, July 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
mep  = 5.44617013e-4;  % electron/proton mass ratio
a0cm = 5.29177249e-9;  % Bohr radius in cm
mp = 1.007276470;      % proton mass in amu
mu = mu/2*mep/(mp*137.0359895);
Q  = Q*1e-24/a0cm^2; 
if (I > 0)
  M(1) = mu*sqrt((I+1)*I)/I;
  if (I > 0.5)
    M(2) = Q/2*sqrt((2*I+3)*(I+1))/(I*(2*I-1));
  else
    M(2) = 0;
  end
else
  M(1) = 0;
  M(2) = 0;
end

M(1) = M(1);
M(2) = M(2);
