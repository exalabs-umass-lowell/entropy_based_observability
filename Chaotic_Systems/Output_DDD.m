function [Ar Cr V] = Output_DDD(Ord)

a = 0.1;
b = 0.1;
c = 14;
k = 0.1;

A = [-k -1 -1 0 k 0 0 0 0;
    1 1 0 0 0 0 0 0 0;
    0 0 -c 0 0 0 0 0 0;
    0 k 0 -k*2 -1 -1  0 k 0;
    0 0 0 1 1 0 0 0 0;
    0 0 0 0 0 -c 0 0 0;
    0 0 0 0 k 0 -k -1 -1;
    0 0 0 0 0 0 1 1 0;
    0 0 0 0 0 0 0 0 -c];



% W matrix for 2nd order terms

W = zeros(9,9,9);
W(1,3,3) = 1/2;
W(3,1,3) = 1/2;
W(4,6,6) = 1/2;
W(6,4,6) = 1/2;
W(7,9,9) = 1/2;
W(9,7,9) = 1/2;


N = 9;
b = ones(N,1);
%%
C = zeros(1,N);
C(1) = 1;

[V,Q] = Arnoldi(inv(A),b,Ord-1);
H = V'*inv(A)*V;
% By using Arnoldi iterations, we got the H matrix and V matrix as follows
Ar = inv(H);
Cr = C*V;