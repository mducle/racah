function Hss = judd_Hss(n,l,M,cellflag)
% Calculates the spin-spin Hamiltonian for the l^n configuration by a chain-calculation

if ~isscalar(n) | ~isscalar(l) | ~isnumeric(n) | ~isnumeric(l) 
  error('n, l must be numerical scalars');
elseif ~isvector(M) | length(M)~=3
  error('M must be a three-component vector');
end


%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
mx = [  1    8/sqrt(3)              0           0           0           0           0;
        0 4*sqrt(14)/3 8*sqrt(11/2)/3           0           0           0           0;
        0            0  4*sqrt(143)/3           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0];
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
m0 = [-12            3              0           0           0           0           0;
        0           -1              2           0           0           0           0;
        0            0              1           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0].*mx;
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
m2 = [-24            1              0           0           0           0           0;
        0            8         -23/11           0           0           0           0;
        0            0         -34/11           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0].*mx;
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
m4 = [-300/11  -100/11              0           0           0           0           0;
        0      -200/11       -325/121           0           0           0           0;
        0            0     -1325/1573           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0;
        0            0              0           0           0           0           0].*mx;

m0 = m0+triu(m0,1)'; m2 = m2+triu(m2,1)'; m4 = m4+triu(m4,1)';

%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
mx = [  1            0              0           1 -sqrt(2/15)           0           0;  % 3P
        0   2*sqrt(14)              0           0   sqrt(2/5)    sqrt(11)           0;  % 3F
        0            0     8/sqrt(55)           0           0   sqrt(2/5)    sqrt(26);  % 3H
        0            0              0           0           0           0           0;  % 1S
        0            0              0           0           0           0           0;  % 1D
        0            0              0           0           0           0           0;  % 1G
        0            0              0           0           0           0           0]; % 1I
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
n0 = [-36            0              0           6          27           0           0;  % 3P
        0          -15              0           0          23          -6           0;  % 3F
        0            0           -132           0           0          39          -5;  % 3H
        0            0              0           0           0           0           0;  % 1S
        0            0              0           0           0           0           0;  % 1D
        0            0              0           0           0           0           0;  % 1G
        0            0              0           0           0           0           0].*mx;
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
n2 = [-72            0              0           2          14           0           0;  % 3P
        0           -1              0           0           6       64/33           0;  % 3F
        0            0             23           0           0     -728/33      -30/11;  % 3H
        0            0              0           0           0           0           0;  % 1S
        0            0              0           0           0           0           0;  % 1D
        0            0              0           0           0           0           0;  % 1G
        0            0              0           0           0           0           0].*mx;
%    '3P'         '3F'           '3H'        '1S'        '1D'        '1G'        '1I'
n4 = [-900/11        0              0       10/11      115/11           0           0;  % 3P
        0        10/11              0           0     -195/11   -1240/363           0;  % 3F
        0            0         130/11           0           0   -3175/363   -375/1573;  % 3H
        0            0              0           0           0           0           0;  % 1S
        0            0              0           0           0           0           0;  % 1D
        0            0              0           0           0           0           0;  % 1G
        0            0              0           0           0           0           0].*mx;

n0 = n0+triu(n0,1)'; n2 = n2+triu(n2,1)'; n4 = n4+triu(n4,1)';

oldm22 = {m0 m2 m4}; oldm11 = {n0 n2 n4}; st_old = racah_states(2,l); s = 1/2; tm = 2; tn = 1;

H = (m0+n0).*M(1) + (m2+n2).*M(2) + (m4+n4).*M(3);
if nargin==4
  Hss.f2 = {m0 m2 m4 n0 n2 n4 H};
else
  Hss = H;
end

if n==2
  return;
elseif n<2 || n>14
  error('number of equivalent electron, n, must be between 2 and 14');
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

  for im = 1:3
    m22{im} = zeros(num_states); m11{im} = zeros(num_states);
    %[i22nz,j22nz] = find(oldm22{im}); [i11nz,j11nz] = find(oldm11{im}); 
    %ittnz = [i22nz i11nz]; jttnz = [j22nz j11nz];
    for i = 1:num_states        % Calculates only the upper triangle
      for j = i:num_states      % Lower triangle is related by: <x|H|x'> = (-1)^((L-S)-(L'-S'))<x'|H|x>
        %icfnz = find(cfp{i}); jcfnz = find(cfp{j});

        %[inz,jnz] = find( (cfp{i}'*cfp{j}).*(oldm22{im}+oldm11{im}) );
	   % With these indices - reduced time from 7.2min to 3.4min for f4 for both t=2 and t=1 matrices
	   % Using separate loops e.g.:
	   %  [inz,jnz] = find( (cfp{i}'*cfp{j}).*oldm22{im} ); 
	   % Reduce to 1.9min for f4
        %for p = 1:length(inz);
        %  pi = inz(p); pj = jnz(p);
        %  %for pj = find(cfp{j})
        %    mat_el = cfp{i}(pi)*cfp{j}(pj) * (-1)^(Sb(pi)+Lb(pi)+s+l+S(j)+L(j)) ...
        %             * sqrt((2*S(i)+1)*(2*S(j)+1)*(2*L(i)+1)*(2*L(j)+1));
        %    m22{im}(i,j) = m22{im}(i,j) + ( mat_el * sixj([S(i) tm S(j); Sb(pj) s Sb(pi)]) ...
        %                                           * sixj([L(i) tm L(j); Lb(pj) l Lb(pi)]) * oldm22{im}(pi,pj) );
        %    m11{im}(i,j) = m11{im}(i,j) + ( mat_el * sixj([S(i) tn S(j); Sb(pj) s Sb(pi)]) ...
        %                                           * sixj([L(i) tn L(j); Lb(pj) l Lb(pi)]) * oldm11{im}(pi,pj) );
        %  %end
        %end

        if abs(S(j)-S(i))<=2 && abs(L(j)-L(i))<=2
        cfplogmat = cfp{i}'*cfp{j}; [inz,jnz] = find( cfplogmat.*oldm22{im} );
        for p = 1:length(inz);
          pi = inz(p); pj = jnz(p);
          mat_el = cfp{i}(pi)*cfp{j}(pj) * (-1)^(Sb(pi)+Lb(pi)+s+l+S(j)+L(j)) ...
                   * sqrt((2*S(i)+1)*(2*S(j)+1)*(2*L(i)+1)*(2*L(j)+1));
          m22{im}(i,j) = m22{im}(i,j) + ( mat_el * sixj([S(i) tm S(j); Sb(pj) s Sb(pi)]) ...
                                                 * sixj([L(i) tm L(j); Lb(pj) l Lb(pi)]) * oldm22{im}(pi,pj) );
        end
        end % if /\S=0,1,2, /\L=0,1,2
	   % Separate out t=1 and t=2 matrices and apply triangular conditions on six-j's (StS') and (LtL')
	   % reduces time to 1.3min for f4
      end
    end

    for i = 1:num_states        % Calculates only the upper triangle
      for j = i:num_states      % Lower triangle is related by: <x|H|x'> = (-1)^((L-S)-(L'-S'))<x'|H|x>
        if abs(S(j)-S(i))<=1 && abs(L(j)-L(i))<=1
        cfplogmat = cfp{i}'*cfp{j}; [inz,jnz] = find( cfplogmat.*oldm11{im} );
        for p = 1:length(inz);
          pi = inz(p); pj = jnz(p);
          mat_el = cfp{i}(pi)*cfp{j}(pj) * (-1)^(Sb(pi)+Lb(pi)+s+l+S(j)+L(j)) ...
                   * sqrt((2*S(i)+1)*(2*S(j)+1)*(2*L(i)+1)*(2*L(j)+1));
          m11{im}(i,j) = m11{im}(i,j) + ( mat_el * sixj([S(i) tn S(j); Sb(pj) s Sb(pi)]) ...
                                                 * sixj([L(i) tn L(j); Lb(pj) l Lb(pi)]) * oldm11{im}(pi,pj) );
        end
        end % if /\S=0,1, /\L=0,1
      end
    end
    m22{im} = m22{im}+m22{im}'.*lTriFact;       m11{im} = m11{im}+m11{im}'.*lTriFact;
    m22{im} = ((n_calc)/(n_calc-2)) .* m22{im}; m11{im} = ((n_calc)/(n_calc-2)) .* m11{im};
  end

  time = toc; if time<60; display(sprintf('Elapsed time was %0.2g s',time));
                     else display(sprintf('Elapsed time was %0.2g min',time/60)); end;

  %if exist('matname') && ~exist(matfilename,'file')
  %  display(sprintf('Saving %s',matfilename));
  %  save(matfilename,'matrix');
  %end

  if nargin==4
    H = (m22{1}+m11{1}).*M(1) + (m22{2}+m11{2}).*M(2) + (m22{3}+m11{3}).*M(3);
    Hss.(sprintf('%c%1g',lower(racah_lconv(l)),n_calc)) = {m22{1} m22{2} m22{3} m11{1} m11{2} m11{3} H};
  end
  
  oldm22 = m22; oldm11 = m11;
  st_old = st_now;

end

if nargin~=4
  Hss = (m22{1}+m11{1}).*M(1) + (m22{2}+m11{2}).*M(2) + (m22{3}+m11{3}).*M(3);
end
