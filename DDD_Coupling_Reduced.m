
function f=DDD_Coupling_Reduced(t,M)



a = 0.1;
b = 0.1;
c = 14;
k = 0.1;

% A = [-k -1 -1 0 k 0 0 0 0;
%     1 1 0 0 0 0 0 0 0;
%     0 0 -c 0 0 0 0 0 0;
%     0 k 0 -k*2 -1 -1  0 k 0;
%     0 0 0 1 1 0 0 0 0;
%     0 0 0 0 0 -c 0 0 0;
%     0 0 0 0 k 0 -k -1 -1;
%     0 0 0 0 0 0 1 1 0;
%     0 0 0 0 0 0 0 0 -c];



A = [-k -1 -1 0 k 0 0 0 0;
    1 1 0 0 0 0 0 0 0;
    0 0 -c 0 0 0 0 0 0;
    0 0 0 0 -1 -1 0 0 0;
    0 k 0 1 1-2*k 0 0 k 0;
    0 0 0 0 0 -c 0 0 0;
    0 0 0 0 0 0 0 -1 -1;
    0 0 0 0 0 0 1 1 0;
    0 0 0 0 0 k 0 0 -k-c];



N = 9;
Ord = size(M,1);

b = ones(N,1);
%%
B = zeros(N,1);
B(1) = 1;

W = zeros(9,9,9);
W(1,3,3) = 1/2;
W(3,1,3) = 1/2;
W(4,6,6) = 1/2;
W(6,4,6) = 1/2;
W(7,9,9) = 1/2;
W(9,7,9) = 1/2;


[V,Q] = Arnoldi(inv(A),b,Ord-1);
H = V'*inv(A)*V;
% By using Arnoldi iterations, we got the H matrix and V matrix as follows

Ar = inv(H);
Br = inv(H)*V'*inv(A)*B;


x = M;

WX = zeros(N,1);

for i = 1:N
    WX(i)= x'*V'*W(:,:,i)*V*x;
end

f = inv(H)*x + inv(H)*V'*inv(A)*WX + Br*5;