function cfp = racah_cfp(n,U,v,S,L,Up,vp,Sp,Lp)
% Calculates the coefficient of fractional parentage for l=3, after the method of Racah IV.

if nargin==3 && ischar(U) && ischar(v)
  if ~ischar(n)
    nc = ['f' char(double(n)+48)];
    np = ['f' char(double(n-1)+48)];
  else
    nc = n;
    np = ['f' char(double(n(2)-1))];
    n = double(n(2)-48);
  end
  stc = racah_states(nc,U);  % Child state
  stp = racah_states(np,v);  % Parent state
  U  = stc{4}; v  = stc{3}; S  = stc{1}; L  = stc{2};
  Up = stp{4}; vp = stp{3}; Sp = stp{1}; Lp = stp{2};
end

nn = n; 
if n>7
  n = 14-nn; nc = ['f' char(double(n+1)+48)]; np = ['f' char(double(n)+48)];
end

%{U(1) U(2) v S L Up(1) Up(2) vp Sp Lp}

if ~ischar(L)
  errmsg = sprintf(['You must specify the orbital angular momentum quantum number as S,P,D,F,G, etc. \n' ...
         'This is to specify certain combinations where the |vUSL> basis is not enough to distinguish \n' ...
	 'between states, eg. 2F6 and 2F7, specified by S=1/2 v=5 U=(31) and L = "F" or "Fp"'],'');
  error(errmsg);
end

Ls = L; Lsp = Lp;
L = racah_lconv(L); Lp = racah_lconv(Lp);

W  = racah_vtow(v,S,3);           % Hard-coded for f-electrons (i.e. l=3)
Wp = racah_vtow(vp,Sp,3); 

ulf  = racah_ulf(U,Ls,Up,Lsp);    % The difference in the L and Lp values only matters in this table.
wupf = racah_wupf(W,U,Wp,Up);

l = 3;                            % Hard-coded for f-electrons (i.e. l=3)

if n==v
  if Sp==(S-(1/2))
    cfp = sqrt( (v+2*S+2)*S   / (v*(2*S+1)) ) * ulf * wupf;
  elseif Sp==(S+(1/2))
    cfp = sqrt( (v-2*S)*(S+1) / (v*(2*S+1)) ) * ulf * wupf;
  else
    cfp = 0;
  end
elseif vp==(v-1)
  if Sp==(S-(1/2))
    cfp = sqrt( (4*l+4-n-v)*(v+2*S+2)*S   / (2*n*(2*l+2-v)*(2*S+1)) ) * ulf * wupf;
  elseif Sp==(S+(1/2))
    cfp = sqrt( (4*l+4-n-v)*(v-2*S)*(S+1) / (2*n*(2*l+2-v)*(2*S+1)) ) * ulf * wupf;
  else
    cfp = 0;
  end
elseif vp==(v+1)
  if Sp==(S-(1/2))
    cfp = sqrt( (n-v)*(4*l+6-v+2*S)*S     / (2*n*(2*l+2-v)*(2*S+1)) ) * ulf * wupf;
  elseif Sp==(S+(1/2))
    cfp = sqrt( (n-v)*(4*l+4-v-2*S)*(S+1) / (2*n*(2*l+2-v)*(2*S+1)) ) * ulf * wupf;
  else
    cfp = 0;
  end
else
  cfp = 0;
end

if mod(v,2)==1                    % v is odd
  cfp = cfp * (-1)^(Sp);
else                              % v is even
  cfp = cfp * (-1)^(Sp+(vp-v)/2);
end

if nn>7
  % Formula for cfp of more than half filled subshell in terms of less than half filled subshells is:
  %
  %   4l+1-n         4l+2-n               S+S'+L+L'-l-1/2 {   (n+1)(2S+1)(2L+1)    }   n            n+1
  % (l       aSL |} l       a'S'L') = (-1)                { ---------------------- } (l  a'S'L' |} l    aSL)
  %                                                       { (4l+2-n)(2S'+1)(2L'+1) }
  %
  % Ref: A.F. Starace, Phys. Rev. A v27, p572 (1983)
  cfp = cfp * ( (-1)^(S+Sp+L+Lp-l-.5) * ((n+1)*(2*S+1)*(2*L+1)/(4*l+2-n)/(2*Sp+1)/(2*Lp+1)) );
end
