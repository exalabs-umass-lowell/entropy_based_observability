function K = Neuron_Lie_matrix(h0,X,N)

K = sym('x',[N N]);


K(1,:) = zeros(1,N);
K(1,3) = 1;
K(2,:) = jacobian(K(1,:)*h0,X);
% M = jacobian(h0,X);
% K(2,:) = M(1,:)*h0;


for i = 3:3*N
    K(i,:) = jacobian(K(i-1,:)*h0,X);   
end








