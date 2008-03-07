function chi = racah_chi(U,Up,Lp,L)
% Calculates the matrix element (U|x(L)|U') x==chi from Racah IV tables.

t1 = [0       0     143       0    -130       0      35       0       0       0       0;
      0       0       0       0       0       1       0       0       0       0       0;
      0       0  -39*sqrt(2)  0   4*sqrt(65)  0       0       0       0       0       0;
      0       0     377     455    -561      49       0    -315     245       0       0;
      0       0      13     -65      55     -75       0     133     -75       0       0;
      0       0       0       1       0       0       0       0       0       0       0;
      0 -13*sqrt(11)  0       0       0   sqrt(39)    0       0       0       0       0;
      0       0       0       0 -13*sqrt(5)   0      30       0       0       0       0;
      0       0       0 12*sqrt(195) 8*sqrt(143) 11*sqrt(42) 0 -4*sqrt(17) 0  0       0;
      0     -52       0      38     -52      88      25     -94       0      25       0;
      0       0  3*sqrt(429)  0 -38*sqrt(65)  0  21*sqrt(85)  0       0       0       0;
      0       0 45*sqrt(78)   0  12*sqrt(11) -12*sqrt(546) 0  0 -8*sqrt(665)  0       0;
    260       0     -25       0      94     104    -181       0     -36       0      40];

t2 = [0    11*sqrt(330)               0               0    76*sqrt(143)           -6644                  0;
      0               0     -8*sqrt(78)  -60*sqrt(39/7)               0            4792                  0;
      0               0               0    -312*sqrt(5)    -48*sqrt(39)            4420      336*sqrt(143);
      1               0               0    12*sqrt(715)    -98*sqrt(33)            -902      336*sqrt(143);
      0               0      5*sqrt(65)    2024/sqrt(7)   20*sqrt(1001)           -2684                  0;
      0     11*sqrt(85)               0 31*sqrt(1309/3)   -20*sqrt(374)           -2024     -48*sqrt(6545);
      0    -25*sqrt(77)               0   103*sqrt(5/3)    -44*sqrt(70)            2680     -48*sqrt(6545);
      0               0     10*sqrt(21)               0    -57*sqrt(33)        -12661/5   -3366*sqrt(34)/5;
      0               0               0               0   18*sqrt(1122)         17336/5   -3366*sqrt(34)/5;
      0               0               0 -52*sqrt(323/23) -494*sqrt(19/23)     123506/23 144*sqrt(21318)/23;
      0               0               0 -336*sqrt(66/23) 73*sqrt(1122/23)     -85096/23 144*sqrt(21318)/23;
      0               0               0   -24*sqrt(190)               0           -4712                  0;
      0               0               0               0   -21*sqrt(385)            -473                  0;
      0               0               0               0               0            1672                  0;
      0               0               0               0               0             220                  0];

%t2 = {t2; zeros(size(t2); zeros(size(t2)); zeros(size(t2))}; 
%t2{2}(6,3)  = 336*sqrt(143);      t2{3}(6,3)  = 336*sqrt(143);      t2{4}(6,3)  = -902;
%t2{2}(6,4)  = 336*sqrt(143);      t2{3}(6,4)  = 336*sqrt(143);      t2{4}(6,4)  = -902;
%t2{2}(6,6)  = -48*sqrt(6545);     t2{3}(6,6)  = -48*sqrt(6545);     t2{4}(6,6)  = 2680;
%t2{2}(6,7)  = -48*sqrt(6545);     t2{3}(6,7)  = -48*sqrt(6545);     t2{4}(6,7)  = 2680;
%t2{2}(6,8)  = -3366*sqrt(34)/5;   t2{3}(6,8)  = -3366*sqrt(34)/5;   t2{4}(6,8)  = 17336/5;
%t2{2}(6,9)  = -3366*sqrt(34)/5;   t2{3}(6,9)  = -3366*sqrt(34)/5;   t2{4}(6,9)  = 17336/5;
%t2{2}(6,10) = 144*sqrt(21318)/23; t2{3}(6,10) = 144*sqrt(21318)/23; t2{4}(6,10) = -85096/23;
%t2{2}(6,11) = 144*sqrt(21318)/23; t2{3}(6,11) = 144*sqrt(21318)/23; t2{4}(6,11) = -85096/23;

t3 = [1               0               0               0           -1408                 0;
      0               0    -88*sqrt(13)               0             -44                 0;
      0               1               0     90*sqrt(11)            1078                 0;
      0               0 55*sqrt(715/27)  -16*sqrt(1001)        -16720/9  -34*sqrt(2618)/9;
      0               0 7*sqrt(15470/27)   64*sqrt(442)         10942/9  -34*sqrt(2618)/9;
      0               0               0   -72*sqrt(462)            -704                 0;
      0               0 34*sqrt(1045/31) -9*sqrt(21945/31)     -2453/31 60*sqrt(74613)/31;
      0               0 -12*sqrt(1785/31) 756*sqrt(85/31)      36088/31 60*sqrt(74613)/31;
      0               0               0    -84*sqrt(33)            -132                 0;
      0               0               0               0        -4268/31 924*sqrt(1995)/31;
      0               0               0               0        11770/31 924*sqrt(1995)/31;
      0               0               0    -99*sqrt(15)           -1067                 0;
      0               0               0               0             528                 0;
      0               0               0               0              22                 0];
%t3 = {t3; zeros(size(t3); zeros(size(t3)); zeros(size(t3))}; 
%t3{2}(5,4)  = 34*sqrt(2618)/9;   t3{3}(5,4)  = 34*sqrt(2618)/9;   t3{4}(5,4)  = 10942/9;
%t3{2}(5,5)  = 34*sqrt(2618)/9;   t3{3}(5,5)  = 34*sqrt(2618)/9;   t3{4}(5,5)  = 10942/9;
%t3{2}(5,7)  = 60*sqrt(74613)/31; t3{3}(5,7)  = 60*sqrt(74613)/31; t3{4}(5,7)  = 36088/31;
%t3{2}(5,8)  = 60*sqrt(74613)/31; t3{3}(5,8)  = 60*sqrt(74613)/31; t3{4}(5,8)  = 36088/31;
%t3{2}(5,10) = 924*sqrt(1995)/31; t3{3}(5,10) = 924*sqrt(1995)/31; t3{4}(5,10) = 11770/31;
%t3{2}(5,11) = 924*sqrt(1995)/31; t3{3}(5,11) = 924*sqrt(1995)/31; t3{4}(5,11) = 11770/31;


if Up==[3 1]         % Use table t2
  if U==[1 0];     j = 1;
  elseif U==[1 1]; j = 2;
  elseif U==[2 0]; j = 3;
  elseif U==[2 1]; j = 4;
  elseif U==[3 0]; j = 5;
  %elseif U==[3 1]; j = 6;
  %elseif U==[3 1]; if L==Lp; j = 6; elseif strncmp(L,Lp,1); j = 7; end;
  elseif U==[3 1]; if L==Lp; j = 6; elseif L(1)==Lp(1); j = 7; end;
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
    chi = t2(i,j);
  else
    chi = 0;
  end

elseif Up==[4 0]     % Use table t3
  if U==[0 0];     j = 1;
  elseif U==[1 0]; j = 2;
  elseif U==[2 0]; j = 3;
  elseif U==[3 0]; j = 4;
  %elseif U==[4 0]; j = 5;
  %elseif U==[4 0]; if L==Lp; j = 5; elseif strncmp(L,Lp,1); j = 6; end;
  elseif U==[4 0]; if L==Lp; j = 5; elseif L(1)==Lp(1); j = 6; end;
  else j = 0; end;
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
    case 'M';  i = 12;
    case 'N';  i = 13;
    case 'Q';  i = 14;
  otherwise; 
    i = 0;
  end

  if i~=0 & j~=0
    chi = t3(i,j);
  else
    chi = 0;
  end

else                % Use table t1
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

      if U==[2 0] & Up==[2 0]; i = 1;
  elseif U==[1 1] & Up==[2 1]; i = 2;
  elseif U==[2 0] & Up==[2 1]; i = 3;
  elseif U==[2 1] & Up==[2 1]; i = 4;    % (21|chi1|21)
  %elseif U==[2 1] & Up==[2 1]; i = 5;    % (21|chi2|21)
  elseif U==[1 0] & Up==[3 0]; i = 6;
  elseif U==[1 1] & Up==[3 0]; i = 7;
  elseif U==[2 0] & Up==[3 0]; i = 8;
  elseif U==[2 1] & Up==[3 0]; i = 9;
  elseif U==[3 0] & Up==[3 0]; i = 10;
  elseif U==[2 0] & Up==[2 2]; i = 11;
  elseif U==[2 1] & Up==[2 2]; i = 12;
  elseif U==[2 2] & Up==[2 2]; i = 13;
  else; i = 0; end
  
  if i==4 & j~=0
    chi = [t1(4,j) t1(5,j)];
  elseif i~=0 & j~=0
    chi = t1(i,j);
  else
    chi = 0;
  end

end
