function F_k = racah_FtoF_k(F,fl)
% Converts the Slater integral parameters F^k to F_k

if ~exist('fl')
  F_k = F(:)'./[1 225 1089 184041/25];
else
  F_k = F(:)'.*[1 225 1089 184041/25];
end
