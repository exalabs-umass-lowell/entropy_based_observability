% This function is named Multi_Neuron, used to simulate the dynamics of
% neuron interactions

% G matrix is the ajacency matrix, showing the connections between
% different neurons

function f = Multi_Neuron_Sym(x)

a = 1;
b = 3.0;
c = 1.0;
d = 5.0;
s = 4.0;
r = 0.006;
x0 = -1.6;
Epsilon = 0.5;
I_ext = 3;

N = size(x,1)/3;
G = ones(N);

for i = 1:N
    G(i,i) = -N+i;
end

for j = 1:N-1
    for m = 1:j
        G(j+1,m) = 0;
    end
end


state_matrix = vec2mat(x,3);
f = sym('x',[3*N,1]);

for i = 3:3:3*N
    x1 = x(i-2);
    y1 = x(i-1);
    z1 = x(i);
    f(i-2) = (y1-a*x1^3+b*x1^2-z1+I_ext + Epsilon*G(i/3,:)*state_matrix(:,1))/r*i;
    f(i-1) = (c-d*x1^2-y1)/r*i;
    f(i)= (s*(x1-x0)-z1);
end

f(3) = f(3)+10;


   
  