function X = rajnak_Xkkl(n,l,k,kp,lp)

s = racah_states(n,l); X = zeros(length(s),length(s));

for kpp = 2:2:6
  X = X + (2*kpp+1)*sixj([k kp kpp; l l lp]) * judd_Vkkk(n,l,[k kp kpp]);
end
