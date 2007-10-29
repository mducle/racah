function phi = racah_phi(U,Up,Lp,L)
% Looks up the value of Racah's (U|phi(L)|Up) matrix elements from Racah IV

t1 = [   0         -11           0           0           0           3           0           0           0           0           0;
         0           0         -11           0          -4           0           7           0           0           0           0;
         0           0           0           1           0           0           0           0           0           0           0;
         0           0   6*sqrt(2)           0    sqrt(65)           0           0           0           0           0           0;
         0           0         -57          63          55        -105           0         -14          42           0           0;
         0    sqrt(11)           0           0           0    sqrt(39)           0           0           0           0           0;
         0           0           0           0   2*sqrt(5)           0           3           0           0           0           0;
         0           0           0   sqrt(195)  -sqrt(143) -2*sqrt(42)           0 -4*sqrt(17)           0           0           0;
         0          83           0         -72          20         -15          42         -28           0           6           0;
         1           0           0           0           0           0           0           0           0           0           0;
         0           0 3*sqrt(429)           0  4*sqrt(65)           0  3*sqrt(85)           0           0           0           0;
       144           0          69           0        -148          72          39           0         -96           0          56;%];
       % Added following check with Wybourne tables in J. Chem. Phys. v31, p340, 1961
         0           0  -6*sqrt(2)           0   -sqrt(65)           0           0           0           0           0           0];


t2 = [   0       sqrt(330)                0     17*sqrt(143)               209                 0;
         0               0     12*sqrt(273)                0              -200                 0;
         1               0      -36*sqrt(5)     -16*sqrt(39)               624     -80*sqrt(143);
         0               0     -3*sqrt(715)      24*sqrt(33)              -616     -80*sqrt(143);
%         1               0      -36*sqrt(5)     -16*sqrt(39)               624     -80*sqrt(143);
%         0               0     -3*sqrt(715)      24*sqrt(33)     -80*sqrt(143)              -616;
         0               0       11*sqrt(7)     4*sqrt(1001)               836                 0;
         0        sqrt(85)  -2*sqrt(1309/3)      sqrt(187/2)           -1353/2   -5*sqrt(6545)/2;
         0        sqrt(77)    -74*sqrt(5/3)    31*sqrt(35/2)             703/2   -5*sqrt(6545)/2;
%         0        sqrt(85)  -2*sqrt(1309/3)      sqrt(187/2)           -1353/2   -5*sqrt(6545)/2;
%         0        sqrt(77)    -74*sqrt(5/3)    31*sqrt(35/2)   -5*sqrt(6545)/2             703/2;
         0               0                0      30*sqrt(33)           -2662/5    528*sqrt(34)/5;
         0               0                0                0             -88/5    528*sqrt(34)/5;
%         0               0                0      30*sqrt(33)           -2662/5    528*sqrt(34)/5;
%         0               0                0                0   -528*sqrt(34)/5             -88/5;
         0               0 -28*sqrt(323/23)      4*sqrt(437)           6652/23 96*sqrt(21318)/23;
         0               0   42*sqrt(66/23)                0          -5456/23 96*sqrt(21318)/23;
%         0               0 -28*sqrt(323/23)      4*sqrt(437)           6652/23 96*sqrt(21318)/23;
%         0               0   42*sqrt(66/23)                0 96*sqrt(21318)/23          -5456/23;
         0               0     -6*sqrt(190)                0              -464                 0;
         0               0                0     -6*sqrt(385)               814                 0;
         0               0                0                0              -616                 0;
         0               0                0                0               352                 0];

t3 = [                   0                   0        2*sqrt(2145);
               11*sqrt(13)         -6*sqrt(26)          9*sqrt(33);
                         0         3*sqrt(455)                   0;
           -4*sqrt(715/27)    -131*sqrt(11/27)      -4*sqrt(11/27);
            sqrt(15470/27)     17*sqrt(238/27)    -17*sqrt(238/27);
                         0        -12*sqrt(21)         3*sqrt(286);
           7*sqrt(1045/31)                   0     3*sqrt(3553/31);
           3*sqrt(1785/31)                   0      75*sqrt(21/31);
                         0        -2*sqrt(119)                   0;
                         0     22*sqrt(105/31)      4*sqrt(627/31);
                         0     -84*sqrt(19/31)     12*sqrt(385/31);
                         0     -84*sqrt(19/31)          sqrt(2530)];

if Up==[3 1]         % Use table t2
      if U==[1 0]; j = 1;
  elseif U==[1 1]; j = 2;
  elseif U==[2 1]; j = 3;
  elseif U==[3 0]; j = 4;
  %elseif U==[3 1]; if length(L)==1; j = 5; elseif strcmp(L(2),'p'); j = 6; end;
  elseif U==[3 1]; if L==Lp; j = 5; elseif strncmp(L,Lp,1); j = 6; end;
  else j = 0; end;
  switch L
    case 'P';  i = 1;
    case 'D';  i = 2;
    case 'F';  i = 3;
    case 'Fp'; i = 4;
    case 'G';  i = 5;
    case 'H';  i = 6;
    case 'Hp'; i = 7;
    case 'I';  i = 8;
    case 'Ip'; i = 9;
    case 'K';  i = 10;
    case 'Kp'; i = 11;
    case 'L';  i = 12;
    case 'M';  i = 13;
    case 'N';  i = 14;
    case 'O';  i = 15;
  otherwise; 
    i = 0;
  end

  if i~=0 & j~=0
    phi = t2(i,j);
  else
    phi = 0;
  end

elseif Up==[4 0]     % Use table t3
      if U==[2 0]; j = 1;
  elseif U==[2 1]; j = 2;
  elseif U==[2 2]; j = 3;
  else j = 0; end
  switch L
    case 'S';  i = 1;
    case 'D';  i = 2;
    case 'F';  i = 3;
    case 'G';  i = 4;
    case 'Gp'; i = 5;
    case 'H';  i = 6;
    case 'I';  i = 7;
    case 'Ip'; i = 8;
    case 'K';  i = 9;
    case 'L';  i = 10;
    case 'Lp'; i = 11;
    case 'N';  i = 12;
  otherwise; 
    i = 0;
  end

  if i~=0 & j~=0
    phi = t3(i,j);
  else
    phi = 0;
  end

else
  switch L
    case 'S'; j = 1;
    case 'P'; j = 2;
    case 'D'; j = 3;
    case 'F'; j = 4;
    case 'G'; j = 5;
    case 'H'; j = 6;
    case 'I'; j = 7;
    case 'K'; j = 8;
    case 'L'; j = 9;
    case 'M'; j = 10;
    case 'N'; j = 11;
  otherwise
    j = 0;
  end

%      if  U==[1 1] & Up==[1 1]; i = 1;
%  elseif  U==[2 0] & Up==[2 0]; i = 2;
%  elseif (U==[1 0] & Up==[2 1]) | (U==[2 1] & Up==[1 0]); i = 3;
%  elseif (U==[2 0] & Up==[2 1]) | (U==[2 1] & Up==[2 0]); i = 4;
%  elseif  U==[2 1] & Up==[2 1]; i = 5;
%  elseif (U==[1 1] & Up==[3 0]) | (U==[3 0] & Up==[1 1]); i = 6;
%  elseif (U==[2 0] & Up==[3 0]) | (U==[3 0] & Up==[2 0]); i = 7;
%  elseif (U==[2 1] & Up==[3 0]) | (U==[3 0] & Up==[2 1]); i = 8;
%  elseif  U==[3 0] & Up==[3 0]; i = 9;
%  elseif (U==[0 0] & Up==[2 2]) | (U==[2 2] & Up==[0 0]); i = 10;
%  elseif (U==[2 0] & Up==[2 2]) | (U==[2 2] & Up==[2 0]); i = 11;
%  elseif  U==[2 2] & Up==[2 2]; i = 12;
%  else; i = 0; end

      if U==[1 1] & Up==[1 1]; i = 1;
  elseif U==[2 0] & Up==[2 0]; i = 2;
%  elseif U==[1 0] & Up==[2 1]; i = 3;
  elseif (U==[1 0] & Up==[2 1]) | (U==[2 1] & Up==[1 0]); i = 3;
  elseif U==[2 0] & Up==[2 1]; i = 4;
  elseif U==[2 1] & Up==[2 1]; i = 5;
  elseif U==[1 1] & Up==[3 0]; i = 6;
  elseif U==[2 0] & Up==[3 0]; i = 7;
  elseif U==[2 1] & Up==[3 0]; i = 8;
  elseif U==[3 0] & Up==[3 0]; i = 9;
  elseif U==[0 0] & Up==[2 2]; i = 10;
  elseif U==[2 0] & Up==[2 2]; i = 11;
  elseif U==[2 2] & Up==[2 2]; i = 12;
  elseif U==[2 1] & Up==[2 0]; i = 13;
  else; i = 0; end
  
  if i~=0 & j~=0
    phi = t1(i,j);
  else
    phi = 0;
  end

end
