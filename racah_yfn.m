function y = racah_yfn(n,v,S,U,vp,Up)
% Looks up the value of Racah function y(fn,v,S,U,v',S,U') from tables in Racah IV

t1 = [0                   0         -6*sqrt(22);
      2                   0                   0;
      0                10/7        2*sqrt(66)/7;
      0        2*sqrt(66)/7                 2/7];

t2 = [              0              0               0  -12*sqrt(33/5)             0              0               0               0;
                    0            6/5               0               0             6              0               0               0;
                    0              0               0   8*sqrt(11/15)             0              0               0               0;
                    0          29/15               0               0          -1/3              0               0               0;
                    0              0             6/7 -8*sqrt(11/147)     4/sqrt(3)              0               0               0;
        8*sqrt(11/15)              0 -8*sqrt(11/147)           -2/21          -4/3              0               0               0;
                    0           -1/3       4/sqrt(3)            -4/3           1/3              0               0               0;
                    0              0               0               0             0              0               0    -12*sqrt(22);
                    0              0               0               0             0  3*sqrt(3/175)  -4*sqrt(33/35)      -sqrt(3/5);
                    0              0               0               0             0        221/140  8*sqrt(11/245)     -sqrt(7/80);
                    0              0               0               0             0 8*sqrt(11/245)             2/7               0;
                    0              0               0               0             0    -sqrt(7/80)               0             1/4];

t3 = [              0              0               0               0             0;
                    0              0               0      9*sqrt(11)             0;
                    0              0       3/sqrt(7)      sqrt(33/7)   -2*sqrt(21);
                    0              0               0     -sqrt(55/3)             0;
                    0           -1/3               0               0          -5/3;
                    0              0             5/7  5*sqrt(11/147)     2/sqrt(3);
          -sqrt(55/3)              0  5*sqrt(11/147)           -4/21          -2/3;
                    0           -5/3       2/sqrt(3)            -2/3          -1/3];

t4 = [              0              0               0      36/sqrt(5)             0    -36*sqrt(2);
                    0      3/sqrt(2)               0               0   3*sqrt(5)/2    -sqrt(39/8);
                    0              0             3/7   -11*sqrt(6)/7    -4*sqrt(3)              0;
        3*sqrt(33/10)              0   3*sqrt(33/98)    3/7/sqrt(11)  -3/2/sqrt(2)   3/2/sqrt(22);
                    0              0               0     43/sqrt(30)             0      4*sqrt(3);
                    0           -5/6               0               0 -5*sqrt(5/72)   -sqrt(13/48);
                    0              0            11/7   -11/7/sqrt(6)     4/sqrt(3)              0;
          43/sqrt(30)              0    11/7/sqrt(6)          25/231 29/6/sqrt(22)   1/22/sqrt(2);
                    0  -5*sqrt(5/72)       4/sqrt(3)   29/6/sqrt(22)         -1/12   1/4/sqrt(11);
            4*sqrt(3)   -sqrt(13/48)               0    1/22/sqrt(2)  1/4/sqrt(11)           1/44];

t5 = [              0              0               0;
                    0              0     -6*sqrt(11);
                    0   -2*sqrt(2/7)    2*sqrt(33/7)];

t6 = [              0              0               0   -48*sqrt(2/5)             0            -36;
                    0      sqrt(6/5)               0               0       sqrt(3)  3*sqrt(13/10);
                    0              0               0     46*sqrt(15)             0     -8*sqrt(6);
                    0   11/3/sqrt(5)               0               0 -19/3/sqrt(2)    sqrt(13/60);
                    0              0    -6*sqrt(2)/7   -22/7/sqrt(3)   8*sqrt(2/3)              0;
         -sqrt(110/3)              0    sqrt(22/147) -16/21/sqrt(11)   5/3/sqrt(2)     1/sqrt(22);
                    0     -sqrt(5)/3     4*sqrt(2/3)    4/3/sqrt(11)   1/3/sqrt(2)    -1/sqrt(22)];

t7 = [              0              0      6/sqrt(55)    2*sqrt(42/5)  6*sqrt(2/55);
                    0              0   -61/sqrt(770)     8*sqrt(3/5)  -6/sqrt(385);
                    0     3*sqrt(22)       sqrt(2/7)        -sqrt(3)     1/sqrt(7);
        -4*sqrt(33/5)              0     -1/sqrt(22)               0    2/sqrt(11)];

t8 = [              0              0               0      2*sqrt(10)             0;
                    0              0    -16/sqrt(77)      -2*sqrt(6)  6*sqrt(2/77);
                    0      -sqrt(66)       sqrt(6/7)               1     sqrt(3/7)];

t9 = [              0              0    -12*sqrt(11);
                    0     6*sqrt(33)               0;
           -sqrt(5/7)   2*sqrt(11/7)              -1];

switch n
  case 2;
    if S==0 && vp==2 && v==2
      if U==[2 0] & Up==[2 0]
        y = 2;
      else
        y = 0;
      end
    else
      y = 0;
    end
  case 3;
    if S==1/2 && vp==3
          if Up==[1 1]; j = 1;
      elseif Up==[2 0]; j = 2;
      elseif Up==[2 1]; j = 3;
      else; j = 0; end
          if v==1 & U==[1 0]; i = 1;
      elseif v==3 & U==[1 1]; i = 2;
      elseif v==3 & U==[2 0]; i = 3;
      elseif v==3 & U==[2 1]; i = 4;
      else; i = 0; end
      
      if i~=0 & j~=0
        y = t1(i,j);
      else
        y = 0;
      end
    else
      y = 0;
    end
  case 4;
    if vp==4
      if S==1
            if Up==[1 0]; j = 1;
        elseif Up==[1 1]; j = 2;
        elseif Up==[2 0]; j = 3;
        elseif Up==[2 1]; j = 4;
        elseif Up==[3 0]; j = 5;
        else; j = 0; end
	    if v==2 & U==[1 0]; i = 1;
	elseif v==2 & U==[1 1]; i = 2;
	elseif v==4 & U==[1 0]; i = 3;
	elseif v==4 & U==[1 1]; i = 4;
	elseif v==4 & U==[2 0]; i = 5;
	elseif v==4 & U==[2 1]; i = 6;
	elseif v==4 & U==[3 0]; i = 7;
        else; i = 0; end
      elseif S==0
            if Up==[2 0]; j = 6;
        elseif Up==[2 1]; j = 7;
        elseif Up==[2 2]; j = 8;
        else; j = 0; end
            if v==0 & U==[0 0]; i = 8;
        elseif v==2 & U==[2 0]; i = 9;
        elseif v==4 & U==[2 0]; i = 10;
        elseif v==4 & U==[2 1]; i = 11;
        elseif v==4 & U==[2 2]; i = 12;
        else; i = 0; end
      else
        i = 0; j = 0;
      end
      if i~=0 & j~=0
        y = t2(i,j);
      else
        y = 0;
      end
    else
      y = 0;
    end
  case 5;
    if S==3/2 && vp==5
          if Up==[1 0]; j = 1;
      elseif Up==[1 1]; j = 2;
      elseif Up==[2 0]; j = 3;
      elseif Up==[2 1]; j = 4;
      elseif Up==[3 0]; j = 5;
      else; j = 0; end
          if v==3 & U==[0 0]; i = 1;
      elseif v==3 & U==[1 0]; i = 2;
      elseif v==3 & U==[2 0]; i = 3;
      elseif v==5 & U==[1 0]; i = 4;
      elseif v==5 & U==[1 1]; i = 5;
      elseif v==5 & U==[2 0]; i = 6;
      elseif v==5 & U==[2 1]; i = 7;
      elseif v==5 & U==[3 0]; i = 8;
      else; i = 0; end

      if i~=0 & j~=0
        y = t3(i,j);
      else
        y = 0;
      end
    elseif S==1/2 && vp==5
          if Up==[1 0]; j = 1;
      elseif Up==[1 1]; j = 2;
      elseif Up==[2 0]; j = 3;
      elseif Up==[2 1]; j = 4;
      elseif Up==[3 0]; j = 5;
      elseif Up==[3 1]; j = 6;
      else; j = 0; end
          if v==1 & U==[1 0]; i = 1;
      elseif v==3 & U==[1 1]; i = 2;
      elseif v==3 & U==[2 0]; i = 3;
      elseif v==3 & U==[2 1]; i = 4;
      elseif v==5 & U==[1 0]; i = 5;
      elseif v==5 & U==[1 1]; i = 6;
      elseif v==5 & U==[2 0]; i = 7;
      elseif v==5 & U==[2 1]; i = 8;
      elseif v==5 & U==[3 0]; i = 9;
      elseif v==5 & U==[3 1]; i = 10;
      else; i = 0; end

      if i~=0 & j~=0
        y = t4(i,j);
      else
        y = 0;
      end
    else
      y = 0;
    end
  case 6;
    if S==2 && v==4 && vp==6
          if Up==[1 1]; j = 1;
      elseif Up==[2 0]; j = 2;
      elseif Up==[2 1]; j = 3;
      else; j = 0; end
          if U==[0 0]; i = 1;
      elseif U==[1 0]; i = 2;
      elseif U==[2 0]; i = 3;
      else; i = 0; end
      if i~=0 & j~=0
        y = t5(i,j);
      else
        y = 0;
      end
    elseif S==1 && vp==6
          if Up==[1 0]; j = 1;
      elseif Up==[1 1]; j = 2;
      elseif Up==[2 0]; j = 3;
      elseif Up==[2 1]; j = 4;
      elseif Up==[3 0]; j = 5;
      elseif Up==[3 1]; j = 6;
      else; j = 0; end
          if v==2 & U==[1 0]; i = 1;
      elseif v==2 & U==[1 1]; i = 2;
      elseif v==4 & U==[1 0]; i = 3;
      elseif v==4 & U==[1 1]; i = 4;
      elseif v==4 & U==[2 0]; i = 5;
      elseif v==4 & U==[2 1]; i = 6;
      elseif v==4 & U==[3 0]; i = 7;
      else; i = 0; end
      if i~=0 & j~=0
        y = t6(i,j);
      else
        y = 0;
      end
    elseif S==0 && vp==6
          if Up==[0 0]; j = 1;
      elseif Up==[1 0]; j = 2;
      elseif Up==[2 0]; j = 3;
      elseif Up==[3 0]; j = 4;
      elseif Up==[4 0]; j = 5;
      else; j = 0; end
          if v==2 & U==[2 0]; i = 1;
      elseif v==4 & U==[2 0]; i = 2;
      elseif v==4 & U==[2 1]; i = 3;
      elseif v==4 & U==[2 2]; i = 4;
      else; i = 0; end
      if i~=0 & j~=0
        y = t7(i,j);
      else
        y = 0;
      end
    else
      y = 0;
    end
  case 7;
    if S==1/2 && v==3 && vp==7
          if Up==[0 0]; j = 1;
      elseif Up==[1 0]; j = 2;
      elseif Up==[2 0]; j = 3;
      elseif Up==[3 0]; j = 4;
      elseif Up==[4 0]; j = 5;
      else; j = 0; end
          if U==[1 1]; i = 1;
      elseif U==[2 0]; i = 2;
      elseif U==[2 1]; i = 3;
      else; i = 0; end
      if i~=0 & j~=0
        y = t8(i,j);
      else
        y = 0;
      end
    elseif S==3/2 && v==3 && vp==7
          if Up==[2 0]; j = 1;
      elseif Up==[2 1]; j = 2;
      elseif Up==[2 2]; j = 3;
      else; j = 0; end
          if U==[0 0]; i = 1;
      elseif U==[1 0]; i = 2;
      elseif U==[2 0]; i = 3;
      else; i = 0; end
      if i~=0 & j~=0
        y = t9(i,j);
      else
        y = 0;
      end
    else
      y = 0;
    end

  otherwise
    y = 0;
end
