function T11 = judd_T11(n,l,nold,oldz)
% Calculates the z_i operators (i=5-11,c) that forms the T^(11) tensor representing the spin-other-orbit Hamiltonian 

if ~isscalar(n) | ~isscalar(l) | ~isnumeric(n) | ~isnumeric(l) 
  error('n, l must be numerical scalars');
end

%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
zx = [  1          0           0         1 1/sqrt(1080)             0               0;  % 3P
        0   sqrt(14)           0         0  sqrt(2/405)   1/sqrt(891)               0;  % 3F
        0          0  1/sqrt(55)         0            0 1/sqrt(98010) 1/sqrt(1019304);  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0]; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z5 = [  1          0           0         0          210             0               0;  % 3P
        0         -4           0         0          120          -264               0;  % 3F
        0          0          55         0            0         -2310            6006;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z6 = [  1          0           0         0           -6             0               0;  % 3P
        0          0           0         0           -6           -12               0;  % 3F
        0          0          -1         0            0          -186            -546;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z7 = [  1          0           0         0         -105             0               0;  % 3P
        0         -4           0         0          -60           132               0;  % 3F
        0          0          55         0            0          1155           -3003;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z8 = [  1          0           0         0            3             0               0;  % 3P
        0          0           0         0            3             6               0;  % 3F
        0          0          -1         0            0            93             273;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z9 = [  1          0           0         0           70             0               0;  % 3P
        0          1           0         0          -35            77               0;  % 3F
        0          0          55         0            0          -770            2002;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z10= [  0          0           0         0            6             0               0;  % 3P
        0          0           0         0          -12           -24               0;  % 3F
        0          0           0         0            0           186             546;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z11= [  0          0           0         0           12             0               0;  % 3P
        0          0           0         0            0             0               0;  % 3F
        0          0           0         0            0            36            -252;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
%zc = (13/3)z12 - 40z13 + (4/3)z14
%zc = [  1          0           0         2          -90             0               0;  % 3P
%        0          1           0         0           45           -99               0;  % 3F
%        0          0          55         0            0           990           -2574;  % 3H
%        0          0           0         0            0             0               0;  % 1S
%        0          0           0         0            0             0               0;  % 1D
%        0          0           0         0            0             0               0;  % 1G
%        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z12= [  1          0           0        16          -90*(-1)        0               0;  % 3P
        0         -1           0         0           45*(-1)      -99*(-1)          0;  % 3F
        0          0          55*(-1)    0            0           990*(-1)       2574;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z13= [  1          0           0         2          -90             0               0;  % 3P
        0          1           0         0           45           -99               0;  % 3F
        0          0          55         0            0           990           -2574;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z14= [ -2          0           0        10         -180*(-1)        0               0;  % 3P
        0         -2           0         0           90*(-1)     -198*(-1)          0;  % 3F
        0          0         110*(-1)    0            0          1980*(-1)       2574*2; % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z12s=[  8          0           0        16          -90*8           0               0;  % 3P
        0          8           0         0           45*8         -99*8             0;  % 3F
        0          0          55*8       0            0           990*8         -2574*8;% 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I
%    '3P'       '3F'        '3H'      '1S'         '1D'          '1G'            '1I'
z14s=[ 10          0           0        20         -180*5           0               0;  % 3P
        0         10           0         0           90*5        -198*5             0;  % 3F
        0          0         110*5       0            0          1980*5         -2574*10;% 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I

z5 = z5+triu(z5,1)';  z6 = z6+triu(z6,1)';    z7 = z7+triu(z7,1)';   z8 = z8+triu(z8,1)';
z9 = z9+triu(z9,1)'; z10 = z10+triu(z10,1)'; z11 = z11+triu(z11,1)'; %zc = zc+triu(zc,1)';
z12= z12+triu(z12,1)'; z13 = z13+triu(z13,1)'; z14 = z14+triu(z14,1)'; %zc = zc+triu(zc,1)';
z12s=z12s+triu(z12s,1)'; z14s = z14s+triu(z14s,1)'; % One-particle part of the operators.

if ~exist('nold') || ~exist('oldz')
  nold = 2;
  %oldz = {z5 z6 z7 z8 z9 z10 z11 zc};
  %T11.f2 = {z5 z6 z7 z8 z9 z10 z11 zc};
  oldz = {z5 z6 z7 z8 z9 z10 z11 z12 []  z14;... %};   % All the 2-particle parts of the operators
          [] [] [] [] []  []  [] z12s z13 z14s}; % All the 1-particle parts of the operators
          %[] [] [] [] []  []  [] z12s z13 z14s}; % All the 1-particle parts of the operators
  %T11.f2 = {z5 z6 z7 z8 z9 z10 z11 z12 z13 z14};
  T11.f2 = {z5 z6 z7 z8 z9 z10 z11 z12+z12s z13 z14+z14s};
elseif isstruct(oldz)
  T11 = oldz;
  oldz = T11.(sprintf('%c%1g',lower(racah_lconv(l)),nold));
end

st_old = racah_states(nold,l); s = 1/2; t = 1;

if n==2
  return;
elseif n>(2*(2*l+1))
  error('number of equivalent electrons, n>2(2l+1), greater than allowed orbital angular momenta');
elseif n<2
  error('number of equivalent electrons is invalid (1, 0 or negative!)');
end

for n_calc = (nold+1):n
    
  display(sprintf('Calculating matrix for %c^%1g',lower(racah_lconv(l)),n_calc)); tic

  st_now = racah_states(n_calc,l); num_states = length(st_now); clear cfp;
  for i = 1:num_states
    S(i) = st_now{i}{1}; Ls = st_now{i}{2}; v  = st_now{i}{3}; U  = st_now{i}{4};
    L(i) = racah_lconv(Ls); LSdiff(i) = L(i)-S(i);
    [coefs,parents,ind_parents]= racah_parents(n_calc,l,v,U,S(i),Ls);
    cfp{i} = zeros(1,length(st_old)); cfp{i}(ind_parents) = coefs;
  end
  clear lTriFact;
  for i = 1:num_states
    lTriFact(:,i) = LSdiff(i)-LSdiff;
  end
  lTriFact = tril((-1).^lTriFact,-1);
  num_old = length(st_old);
  for i = 1:num_old
    Sb(i) = st_old{i}{1}; Lb(i) = racah_lconv(st_old{i}{2});
  end  

  %zm = zeros(num_states,num_states,num_old,num_old);
  for im = [1:8 10] %1:length(oldz)
    z{1,im} = sparse(num_states,num_states);
    for i = 1:num_states        % Calculates only the upper triangle
      for j = i:num_states      % Lower triangle is related by: <x|H|x'> = (-1)^((L-S)-(L'-S'))<x'|H|x>
        if abs(S(j)-S(i))<=1 && abs(L(j)-L(i))<=1
	  if im==1; zm{i,j} = sparse(num_old,num_old); end
          [inz,jnz] = find( cfp{i}'*cfp{j}.*oldz{1,im} );
          for p = 1:length(inz);
            pi = inz(p); pj = jnz(p);
	    if zm{i,j}(pi,pj)==0
	      zm{i,j}(pi,pj) = ( cfp{i}(pi)*cfp{j}(pj) * (-1)^(Sb(pi)+Lb(pi)+s+l+S(j)+L(j)) ...
                                * sqrt((2*S(i)+1)*(2*S(j)+1)*(2*L(i)+1)*(2*L(j)+1)) ...
                                * sixj([S(i) t S(j); Sb(pj) s Sb(pi)]) * sixj([L(i) t L(j); Lb(pj) l Lb(pi)]) );
            end
	    z{1,im}(i,j) = z{1,im}(i,j) + zm{i,j}(pi,pj)*oldz{1,im}(pi,pj);
	    %z{im}(i,j) = z{im}(i,j) + ( cfp{i}(pi)*cfp{j}(pj) * (-1)^(Sb(pi)+Lb(pi)+s+l+S(j)+L(j)) ...
            %                            * sqrt((2*S(i)+1)*(2*S(j)+1)*(2*L(i)+1)*(2*L(j)+1)) ...
            %                            * sixj([S(i) t S(j); Sb(pj) s Sb(pi)]) ...
            %                            * sixj([L(i) t L(j); Lb(pj) l Lb(pi)]) * oldz{im}(pi,pj) );
          end
        end % if /\S=0,1,2, /\L=0,1,2
      end
    end
    z{1,im} = z{1,im}+z{1,im}'.*lTriFact;
    z{1,im} = ((n_calc)/(n_calc-2)) .* z{1,im}; 
  end

  for im = 8:10
    z{2,im} = racah_chaincal(n_calc,l,n_calc-1,oldz{2,im});
  end

  time = toc; if time<60; display(sprintf('Elapsed time was %0.2g s',time));
                     else display(sprintf('Elapsed time was %0.2g min',time/60)); end;

  %if exist('matname') && ~exist(matfilename,'file')
  %  display(sprintf('Saving %s',matfilename));
  %  save(matfilename,'matrix');
  %end

  oldz = z;
  st_old = st_now;

  T11.(sprintf('%c%1g',lower(racah_lconv(l)),n_calc)) = {z{1,1:7} z{1,8}+z{2,8} z{2,9} z{1,10}+z{2,10}};

end

% Time to to calculate all 8 matrices: f5=2.3min f6=8.4min f7=11min
