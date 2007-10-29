function x = racah_xwu(W,U,Up)
% Looks the Racah function x(W,U,Up) from tables from Racah IV

t1 = {            0               0    12*sqrt(455);
                  0            -6/7    6*sqrt(66)/7;
       12*sqrt(455)    6*sqrt(66)/7          [3/7 0]};  % t1(3,3) = 0 for chi2


t2 = {            0               0               0               0   -20*sqrt(143);
                  0               0               0    10*sqrt(182)              10;
                  0               0            -8/7    4*sqrt(33)/7       4*sqrt(3);
                  0    10*sqrt(182)    4*sqrt(33)/7         [4/7 3]               2; % t2(4,4) = 3 for chi2
      -20*sqrt(143)              10       4*sqrt(3)               2               2};

t3 = {        3/14   3*sqrt(55)/7  -3*sqrt(5/28);
      3*sqrt(55)/7       -[6/7 3]      3/sqrt(7);   % t3(2,2) = -3 for chi2
     -3*sqrt(5/28)      3/sqrt(7)            3/2};

t4 = {            0               0               0               0      5*sqrt(143)   -15*sqrt(429);
                  0               0               0 14*sqrt(910/11)       2*sqrt(10)   2*sqrt(39)/11;
                  0               0             2/7   -10*sqrt(6)/7          sqrt(3)     9*sqrt(3/7);
                  0 14*sqrt(910/11)   -10*sqrt(6)/7    [-1/7 12/11]     5*sqrt(2/11)    3*sqrt(2)/11;  % t4(4,4) = 12/11 for chi2
        5*sqrt(143)      2*sqrt(10)         sqrt(3)    5*sqrt(2/11)             -1/2    3/2/sqrt(11);
      -15*sqrt(429)   2*sqrt(39)/11     9*sqrt(3/7)    3*sqrt(2)/11     3/2/sqrt(11)            1/22};

t5 = {            0               0               0               0   -30*sqrt(143);
                  0               0               0   -3*sqrt(1430)    9*sqrt(1430);
                  0               0            6/11  -3*sqrt(42/11)    9*sqrt(2)/11;
                  0   -3*sqrt(1430)  -3*sqrt(42/11)              -3      1/sqrt(11);
      -30*sqrt(143)    9*sqrt(1430)    9*sqrt(2)/11      1/sqrt(11)            3/11};

if W==[2 0 0]
  if U==[2 0] & Up==[2 0]
    x = 2;
  else
    x = 0;
  end

elseif W==[1 1 1]
  if U==[2 0] & Up==[2 0]
    x = 0;
  else
    x = 0;
  end

elseif W==[2 1 0]

      if U==[1 1]; i = 1;
  elseif U==[2 0]; i = 2;
  elseif U==[2 1]; i = 3;
  else; i = 0; end

      if Up==[1 1]; j = 1;
  elseif Up==[2 0]; j = 2;
  elseif Up==[2 1]; j = 3;
  else; i = 0; end

  if i~=0 & j~=0
    x = t1{i,j};
  else
    x = 0;
  end

elseif W==[2 1 1]

      if U==[1 0]; i = 1;
  elseif U==[1 1]; i = 2;
  elseif U==[2 0]; i = 3;
  elseif U==[2 1]; i = 4;
  elseif U==[3 0]; i = 5;
  else; i = 0; end

      if Up==[1 0]; j = 1;
  elseif Up==[1 1]; j = 2;
  elseif Up==[2 0]; j = 3;
  elseif Up==[2 1]; j = 4;
  elseif Up==[3 0]; j = 5;
  else; j = 0; end

  if i~=0 & j~=0
    x = t2{i,j};
  else
    x = 0;
  end

elseif W==[2 2 0]

      if U==[2 0]; i = 1;
  elseif U==[2 1]; i = 2;
  elseif U==[2 2]; i = 3;
  else; i = 0; end

      if Up==[2 0]; j = 1;
  elseif Up==[2 1]; j = 2;
  elseif Up==[2 2]; j = 3;
  else; i = 0; end

  if i~=0 & j~=0
    x = t3{i,j};
  else
    x = 0;
  end

elseif W==[2 2 1]

      if U==[1 0]; i = 1;
  elseif U==[1 1]; i = 2;
  elseif U==[2 0]; i = 3;
  elseif U==[2 1]; i = 4;
  elseif U==[3 0]; i = 5;
  elseif U==[3 1]; i = 6;
  else; i = 0; end

      if Up==[1 0]; j = 1;
  elseif Up==[1 1]; j = 2;
  elseif Up==[2 0]; j = 3;
  elseif Up==[2 1]; j = 4;
  elseif Up==[3 0]; j = 5;
  elseif Up==[3 1]; j = 6;
  else; j = 0; end

  if i~=0 & j~=0
    x = t4{i,j};
  else
    x = 0;
  end


elseif W==[2 2 2]

      if U==[0 0]; i = 1;
  elseif U==[1 0]; i = 2;
  elseif U==[2 0]; i = 3;
  elseif U==[3 0]; i = 4;
  elseif U==[4 0]; i = 5;
  else; i = 0; end

      if Up==[0 0]; j = 1;
  elseif Up==[1 0]; j = 2;
  elseif Up==[2 0]; j = 3;
  elseif Up==[3 0]; j = 4;
  elseif Up==[4 0]; j = 5;
  else; j = 0; end

  if i~=0 & j~=0
    x = t5{i,j};
  else
    x = 0;
  end

else 
  x = 0;
end
