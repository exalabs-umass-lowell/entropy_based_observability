% Order reduction of Rossler system with Quadratic method
% Parameters: a = 0.25, b = 0.2, c = 1
% W matrix is a 3 by 3 by 3 matrix, with only x(1,3,3) = x(3,1,3) = 1/2,
% other elements are all zero.

function y = Quadratic_Rossler(t,x)

% parameters
A = [0,-1,-1;1,0.250,0;0,0,-1];
b = ones(3,1);
B = [1;0;0];
C = [1 0 0];




[V,Q] = Arnoldi(inv(A),b,1);
H = V'*inv(A)*V;

Ar = inv(H);
Br = inv(H)*V'*inv(A)*B;
Cr = C*V;



W = zeros(3,3,3);
W(3,1,3) = 1/2;
W(1,3,3) = 1/2;

% Wr matrix is the reduced W matrix
Wr = zeros(3,1);
Wr(1) = x'*V'*W(:,:,1)*V*x;
Wr(2) = x'*V'*W(:,:,2)*V*x;
Wr(3) = x'*V'*W(:,:,3)*V*x;

y = inv(H)*x(:) + inv(H)*V'*inv(A)*Wr + Br*sin(t);


