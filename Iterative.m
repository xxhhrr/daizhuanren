function [E,FAI]= diedai(sol,m,n,L,N,w)  
%m,n为当前元素的坐标
fai = sol(m,n);
FAI = sol(m,n);
Q0 = cost(sol,L,N,w);
E0 = Q0(1) * w(1) +  Q0(2) * w(2) + Q0(3) * w(3) + Q0(4) * w(4);
E1 = 0;
E2 = 0;
E3 = 0;
E4 = 0;
E = E0;
if fai ~= 0
    sol(m,n) = 0;
    Q1 = cost(sol,L,N,w);
    E1 = Q1(1) * w(1) +  Q1(2) * w(2) + Q1(3) * w(3) + Q1(4) * w(4);
    if E1 < E
        E = E1;
        FAI = 0;
    end
end
if fai ~= pi/2
    sol(m,n) = pi/2;
    Q2 = cost(sol,L,N,w);
    E2 = Q2(1) * w(1) +  Q2(2) * w(2) + Q2(3) * w(3) + Q2(4) * w(4);
    if E2 < E
        E = E2;
        FAI = pi/2;
    end
end
if fai ~= pi
    sol(m,n) = pi;
    Q3 = cost(sol,L,N,w);
    E3 = Q3(1) * w(1) +  Q3(2) * w(2) + Q3(3) * w(3) + Q3(4) * w(4);
    if E3 < E
        E = E3;
        FAI = pi;
    end
end
if fai ~= 3*pi/2
    sol(m,n) = 3*pi/2;
    Q4 = cost(sol,L,N,w);
    E4 = Q4(1) * w(1) +  Q4(2) * w(2) + Q4(3) * w(3) + Q4(4) * w(4);
    if E4 < E
        E = E4;
        FAI = 3*pi/2;
    end
end
