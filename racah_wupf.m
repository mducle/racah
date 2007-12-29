function wupf = racah_wupf(W,U,Wp,Up)
% Looks up the racah expression (WU | W'U'+f) given W,U,W',U'
% Tables from Racah IV, Phys. Rev. 76 (1949) 1352 for tableIIIa and tableIIIb
% and from Wybourne, J. Chem. Phys. 36 (1961) 2295 for table IIa

%tableIIIa = [0     1     0     0     0     0     0     0     0     0     0;
%             1     0    1/3   2/3    1     0     0     0     0     0     0;
%             0     1     0     0     0   3/35   2/5  18/35  2/5   3/5    0;
%             0     1     0     0     0     0  -1/10  -9/10   0   3/35  32/35;
%             0     1     0     0     0     0     0     0   2/15  9/35  64/105;
%             0     0     1     0     0     0    -1     0     0     0     0;
%             0     0    2/3  -1/3    0   -1/7  -3/8  27/56   0     0     0;
%             0     0    2/9  -7/9    0     0    1/8  -7/8    0     0     0;
%             0     0     1     0     1     0     0     0     0     0     0;
%             0     0    7/9   2/9    1     0     0     0     0     0     0;
%             0     0     0     1     1     0     0     0     0     0     0;
%             0     0     0     0     0   27/35 -9/40 1/280 -3/5   2/5    0;
%             0     0     0     0     0     0    9/10 -1/10   0   32/35 -3/35;
%             0     0     0     0     0     0    7/8   1/8   1/3   2/7  -8/21;
%             0     0     0     0     0     0     0     1   3/16 -25/112 33/56;
%             0     0     0     0     0     0     0     1     0    1/7   6/7;
%             0     0     0     0     0     0     0     0   8/15 -16/35  1/105;
%             0     0     0     0     0     0     0     0   -1/8  27/56 11/28;
%             0     0     0     0     0     0     0     0     0     0     1];

% W' (000) (100)    (110)   (200)      (111)             (210)                   (211)                          (220)
% U' (00)  (10)  (10) (11)  (20)  (00) (10) (20)    (11) (20) (21)     (10) (11)  (20) (21) (30)            (20) (21) (22)
t = [0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0            0     0     0;        % 00 000
     1     0    1/3   2/3    1     0     0     0     0     0     0     0     0     0     0     0            0     0     0;        % 10 100
     0     1     0     0     0   3/35   2/5  18/35  2/5   3/5    0     0     0     0     0     0            0     0     0;        % 10 110
     0     1     0     0     0     0  -1/10  -9/10   0   3/35  32/35   0     0     0     0     0            0     0     0;        % 11 
     0     1     0     0     0     0     0     0   2/15  9/35  64/105  0     0     0     0     0            0     0     0;        % 20 200
     0     0     1     0     0     0    -1     0     0     0     0    [1     0     0     0     0]./1        0     0     0;        % 00 111
     0     0    2/3  -1/3    0   -1/7  -3/8  27/56   0     0     0   [-1     8    15     0     0]./24       0     0     0;        % 10
     0     0    2/9  -7/9    0     0    1/8  -7/8    0     0     0    [1   -56   135  2560  3080]./5832     0     0     0;        % 20
     0     0     1     0     1     0     0     0     0     0     0   [-7     0    15    20     0]./42    9/14 -5/14     0;        % 11 210
     0     0    7/9   2/9    1     0     0     0     0     0     0   [98   448   270  -500   385]./1701  -2/7   5/7     0;        % 20
     0     0     0     1     1     0     0     0     0     0     0    [0    -7   -60   220   385]./672     [9   880  2695]./3584; % 21
     0     0     0     0     0   27/35 -9/40 1/280 -3/5   2/5    0   [-5    40   -27     0     0]./72       0     0     0;        % 10 211
     0     0     0     0     0     0    9/10 -1/10   0   32/35 -3/35 [35     0   -27    64     0]./126      0     0     0;        % 11
     0     0     0     0     0     0    7/8   1/8   1/3   2/7  -8/21 [-245 -280 -867  -512   616]./2520     0     0     0;        % 20
     0     0     0     0     0     0     0     1   3/16 -25/112 33/56 [0    35   -27  -176    77]./315      0     0     0;        % 21
     0     0     0     0     0     0     0     1     0    1/7   6/7   [0     0    27    64  -224]./315      0     0     0;        % 30
     0     0     0     0     0     0     0     0   8/15 -16/35  1/105  0     0     0     0     0            0     0     0;        % 20 220
     0     0     0     0     0     0     0     0   -1/8  27/56 11/28   0     0     0     0     0            0     0     0;        % 21
     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0            0     0     0;        % 22
     0     0     0     0     0     0     0     0     0     0     0    [8     1     0     0     0]./9       -1     0     0;        % 10 221
     0     0     0     0     0     0     0     0     0     0     0   [35     0    27    -1     0]./63    5/14  9/14     0;        % 11
     0     0     0     0     0     0     0     0     0     0     0    [0    35   -24     4     0]./63     5/7   2/7     0;        % 20
     0     0     0     0     0     0     0     0     0     0     0    [0   385 -1188   100  -343]./2016 [-165 -3888  1323]./5376; % 21
     0     0     0     0     0     0     0     0     0     0     0    [0     0     0     7     2]./9    -5/14  9/14     0;        % 30
     0     0     0     0     0     0     0     0     0     0     0    [0     0     0     1     2]./3        0  -1/6   5/6];       % 31
%                                                                                                           ^     ^     ^
% These columns not complete! Partially reassembled from tabulated cfp listed in Neilson and Foster, 1963 --|-----|-----|
% Table continues below, reconstructed from tabulated cfp and also from the code of Allison and McNulty CPC 8 (1974) 246
% And also using the reciprocity relations (Judd, Op. Techniques in At. Spectr., p.181, eq 7-36).

% W'                       (221)                            %  U  W
% U' (10)    (11)     (20)     (21)     (30)     (31)       %
s = [0        0        0        0        0        0;        % 00 000
     0        0        0        0        0        0;        % 10 100
     0        0        0        0        0        0;        % 10 110
     0        0        0        0        0        0;        % 11 
     0        0        0        0        0        0;        % 20 200
     0        0        0        0        0        0;        % 00 111
     0        0        0        0        0        0;        % 10
     0        0        0        0        0        0;        % 20
     0        0        0        0        0        0;        % 11 210
     0        0        0        0        0        0;        % 20
     0        0        0        0        0        0;        % 21
     4/9      5/9      0        0        0        0;        % 10 211
     1/36     0        15/28    55/126   0        0;        % 11
     0        1/9     -4/21    -44/63    0        0;        % 20
     0       -1/576    3/224    25/1008  539/1152 63/128;   % 21
     0        0        0       -7/99     1/9      81/99;    % 30
    -28/243   20/243   20/63   -55/1701 -110/243  0;        % 20 220
     0        1/16     3/56    -9/28     11/32   -7/32;     % 21
     0        0        0        9/99     0        90/99;    % 22
     5/36    -1/9      3/4      0        0        0;        % 10 221
    -1/18     0       -27/70    176/315  0        0;        % 11
     7/36    -1/5     -35/300   0       -44/90    0;        % 20
     0        11/90    0       -224/495 -5/36     63/220;   % 21
     0        0       -36/210  -80/693  -242/360  9/220;    % 30
     0        0        0        16/165   1/60    -117/132]; % 31



%     2/3      0.74536  0        0        0        0;        % 10 211
%     1/6      0        0.73193  0.66069  0        0;        % 11
%     0        1/3     -0.43644 -0.83571  0        0;        % 20
%     0       -0.041667 0.11573  0.15749  0.68402  0.70156;  % 21
%     0        0        0       -0.26591  1/3      0.90453;  % 30
%    -0.33945  0.28689  0.56344 -0.17982 -0.67281  0;        % 20 220
%     0        1/4      0.23146 -0.56695  0.5863  -0.46771;  % 21
%     0        0        0        0.30151  0        0.95346;  % 22
%     0.37268 -1/3      0.86603  0        0        0;        % 10 221
%    -0.2357   0       -0.62106  0.74748  0        0;        % 11
%     0.44096 -0.44721 -0.34157  0       -0.69921  0;        % 20
%     0        0.3496   0       -0.6727  -0.37268  0.53513;  % 21
%     0        0       -0.41404 -0.33976 -0.81989  0.20226;  % 30
%     0        0        0        0.3114   0.1291  -0.94147]; % 31
%s = s.^2.*sign(s);
%s = [0        0        0        0        0        0        0;        % 00 000
%     0        0        0        0        0        0        0;        % 10 100
%     0        0        0        0        0        0        0;        % 10 110
%     0        0        0        0        0        0        0;        % 11 
%     0        0        0        0        0        0        0;        % 20 200
%     0        0        0        0        0        0        0;        % 00 111
%     0        0        0        0        0        0        0;        % 10
%     0        0        0        0        0        0        0;        % 20
%     0        0        0        0        0        0        0;        % 11 210
%     0        0        0        0        0        0        0;        % 20
%     0.86715  0        0        0        0        0        0;        % 21
%     0        0.66667  0.74536  0        0        0        0;        % 10 211
%     0        0.16667  0        0.73193  0.66069  0        0;        % 11
%     0        0        0.33333 -0.43644 -0.83571  0        0;        % 20
%     0        0       -0.041667 0.11573  0.15749  0.68402  0.70156;  % 21
%     0        0        0        0       -0.26591  0.33333  0.90453;  % 30
%     0       -0.33945  0.28689  0.56344 -0.17982 -0.67281  0;        % 20 220
%     0        0        0.25     0.23146 -0.56695  0.5863  -0.46771;  % 21
%     0        0        0        0        0.30151  0        0.95346;  % 22
%     0        0.37268 -0.33333  0.86603  0        0        0;        % 10 221
%     0       -0.2357   0       -0.62106  0.74748  0        0;        % 11
%     0        0.44096 -0.44721 -0.34157  0       -0.69921  0;        % 20
%     0.49608  0        0.3496   0       -0.6727  -0.37268  0.53513;  % 21
%     0        0        0       -0.41404 -0.33976 -0.81989  0.20226;  % 30
%     0.91287  0        0        0        0.3114   0.1291  -0.94147]; % 31

t = [t s];

%tableIIIb = [  [1     0     0     0     0] ./ 1;
%              [-1     8    15     0     0] ./ 24;
%               [1   -56   135  2560  3080] ./ 5823;
%              [-7     0    15    20     0] ./ 42;
%              [98   448   270  -500   385] ./ 1701;
%               [0    -7   -60   220   385] ./ 672;
%              [-5    40   -27     0     0] ./ 72;
%              [35     0   -27    64     0] ./ 126;
%            [-245  -280  -867  -512   616] ./ 2520;
%               [0    35   -27  -176    77] ./ 315;
%               [0     0    27    64  -224] ./ 315];

%tableIIa =  [  [8     1     0     0     0] ./ 9;
%              [35     0    27    -1     0] ./ 63;
%               [0    35   -24     2     0] ./ 63;
%               [0   385 -1188    10  -343] ./ 2016;
%               [0     0     0     7     2] ./ 9;
%               [0     0     0     1     2] ./ 3];

if W==[0 0 0]
  if U == [0 0]; i = 1; else i = 0; end
elseif W==[1 0 0]
  if U == [1 0]; i = 2; else i = 0; end
elseif W==[1 1 0]
  if U == [1 0]; i = 3; elseif U == [1 1]; i = 4; else; i = 0; end
elseif W==[2 0 0]
  if U == [2 0]; i = 5; else i = 0; end
elseif W==[1 1 1] 
  if (U(2)==0 & U(1)<3); i = U(1)+6; else; i = 0; end
elseif W==[2 1 0] 
  if U==[1 1];     i = 9; 
  elseif U==[2 0]; i = 10; 
  elseif U==[2 1]; i = 11; 
  else; i = 0; end
elseif W==[2 1 1]
  if U==[1 0];     i = 12; 
  elseif U==[1 1]; i = 13; 
  elseif U==[2 0]; i = 14; 
  elseif U==[2 1]; i = 15; 
  elseif U==[3 0]; i = 16; 
  else; i = 0; end
elseif W==[2 2 0]
  if U==[2 0];     i = 17; 
  elseif U==[2 1]; i = 18; 
  elseif U==[2 2]; i = 19; 
  else; i = 0; end
elseif W==[2 2 1]         % Table IIa from Wybourne
  if U==[1 0];     i = 20; 
  elseif U==[1 1]; i = 21; 
  elseif U==[2 0]; i = 22; 
  elseif U==[2 1]; i = 23; 
  elseif U==[3 0]; i = 24; 
  elseif U==[3 1]; i = 25; 
  else; i = 0; end
else
  i = 0;
end

if Wp==[0 0 0]
  if Up == [0 0]; j = 1; else j = 0; end
elseif Wp==[1 0 0]
  if Up == [1 0]; j = 2; else j = 0; end
elseif Wp==[1 1 0]
  if Up == [1 0]; j = 3; elseif Up == [1 1]; j = 4; else; j = 0; end
elseif Wp==[2 0 0]
  if Up == [2 0]; j = 5; else j = 0; end
elseif Wp==[1 1 1] 
  if (Up(2)==0 & Up(1)<3); j = Up(1)+6; else; j = 0; end
elseif Wp==[2 1 0] 
  if Up==[1 1];     j = 9; 
  elseif Up==[2 0]; j = 10; 
  elseif Up==[2 1]; j = 11; 
  else; j = 0; end
elseif Wp==[2 1 1]
  if Up==[1 0];     j = 12; 
  elseif Up==[1 1]; j = 13; 
  elseif Up==[2 0]; j = 14; 
  elseif Up==[2 1]; j = 15; 
  elseif Up==[3 0]; j = 16; 
  else; j = 0; end
elseif Wp==[2 2 0]
  if Up==[2 0];     j = 17; 
  elseif Up==[2 1]; j = 18; 
  elseif Up==[2 2]; j = 19; 
  else; j = 0; end
elseif Wp==[2 2 1]
  if Up==[1 0];     j = 20; 
  elseif Up==[1 1]; j = 21; 
  elseif Up==[2 0]; j = 22; 
  elseif Up==[2 1]; j = 23; 
  elseif Up==[3 0]; j = 24; 
  elseif Up==[3 1]; j = 25; 
  else; j = 0; end
else
  j = 0;
end

if i~=0 & j~=0
wupf = real(sqrt(t(i,j))) - imag(sqrt(t(i,j)));
else
  wupf = 0;
end
