function T11 = judd_T11(n,l)
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
zc = [  1          0           0         2          -90             0               0;  % 3P
        0          1           0         0           45           -99               0;  % 3F
        0          0          55         0            0           990           -2574;  % 3H
        0          0           0         0            0             0               0;  % 1S
        0          0           0         0            0             0               0;  % 1D
        0          0           0         0            0             0               0;  % 1G
        0          0           0         0            0             0               0].*zx; % 1I

z5 = z5+triu(z5,1)';  z6 = z6+triu(z6,1)';    z7 = z7+triu(z7,1)';   z8 = z8+triu(z8,1)';
z9 = z9+triu(z9,1)'; z10 = z10+triu(z10,1)'; z11 = z11+triu(z11,1)'; zc = zc+triu(zc,1)';

oldz = {z5 z6 z7 z8 z9 z10 z11 zc}; st_old = racah_states(2,l); s = 1/2; t = 1;

T11.f2 = {z5 z6 z7 z8 z9 z10 z11 zc};

if n==2
  return;
elseif n>(2*(2*l+1))
  error('number of equivalent electrons, n>2(2l+1), greater than allowed orbital angular momenta');
elseif n<2
  error('number of equivalent electrons is invalid (1, 0 or negative!)');
end

for n_calc = 3:n
    
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
  for i = 1:length(st_old)
    Sb(i) = st_old{i}{1}; Lb(i) = racah_lconv(st_old{i}{2});
  end  

  for im = 1:8
    z{im} = zeros(num_states);
    for i = 1:num_states        % Calculates only the upper triangle
      for j = i:num_states      % Lower triangle is related by: <x|H|x'> = (-1)^((L-S)-(L'-S'))<x'|H|x>
        if abs(S(j)-S(i))<=2 && abs(L(j)-L(i))<=2
          [inz,jnz] = find( cfp{i}'*cfp{j}.*oldz{im} );
          for p = 1:length(inz);
            pi = inz(p); pj = jnz(p);
	    z{im}(i,j) = z{im}(i,j) + ( cfp{i}(pi)*cfp{j}(pj) * (-1)^(Sb(pi)+Lb(pi)+s+l+S(j)+L(j)) ...
                                        * sqrt((2*S(i)+1)*(2*S(j)+1)*(2*L(i)+1)*(2*L(j)+1)) ...
                                        * sixj([S(i) t S(j); Sb(pj) s Sb(pi)]) ...
                                        * sixj([L(i) t L(j); Lb(pj) l Lb(pi)]) * oldz{im}(pi,pj) );
          end
        end % if /\S=0,1,2, /\L=0,1,2
      end
    end
    z{im} = z{im}+z{im}'.*lTriFact;
    z{im} = ((n_calc)/(n_calc-2)) .* z{im}; 
  end

  time = toc; if time<60; display(sprintf('Elapsed time was %0.2g s',time));
                     else display(sprintf('Elapsed time was %0.2g min',time/60)); end;

  %if exist('matname') && ~exist(matfilename,'file')
  %  display(sprintf('Saving %s',matfilename));
  %  save(matfilename,'matrix');
  %end

  T11.(sprintf('%c%1g',lower(racah_lconv(l)),n_calc)) = z;
  
  oldz = z;
  st_old = st_now;

end
