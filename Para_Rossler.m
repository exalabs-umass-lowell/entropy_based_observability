function [V Cr] = Para_Rossler()

A = [0,-1,-1;1,0.250,0;0,0,-1];
b = ones(3,1);
B = [1;0;0];
C = [1 0 0];



[V,Q] = Arnoldi(inv(A),b,1);
H = V'*inv(A)*V;

Ar = inv(H);
Br = inv(H)*V'*B;
Cr = C*V;