function cfp = cowan_cfp(n,U,v,S,L,Up,vp,Sp,Lp)
% Calculates the coefficient of fractional parentage for l=3, after the method of Racah IV.

if ischar(n)
  n = double(n(2)-48);
end

if nargin==3 && ischar(U) && ischar(v)
  switch n
    case 5; cfp = cowan_cfp_f5('f5',U,v);
    case 6; cfp = cowan_cfp_f6('f6',U,v);
    case 7; cfp = cowan_cfp_f7('f7',U,v);
  otherwise
    error('Sorry only n=5,6,7 for l=3 implemented so far.');
  end
else
  switch n
    case 5; cfp = cowan_cfp_f5('f5',U,v,S,L,Up,vp,Sp,Lp);
    case 6; cfp = cowan_cfp_f6('f6',U,v,S,L,Up,vp,Sp,Lp);
    case 7; cfp = cowan_cfp_f7('f7',U,v,S,L,Up,vp,Sp,Lp);
  otherwise
    error('Sorry only n=5,6,7 for l=3 implemented so far.');
  end
end
