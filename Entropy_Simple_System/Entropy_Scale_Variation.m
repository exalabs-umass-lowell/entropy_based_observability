% This code tries to find the variation of entropy for scale_average under
% different scales

% The total number we gonna use in the simulation
N = 100;

% Produce time series
t = linspace(0,20,500);
alpha = linspace(1,4,50);
dim = linspace(0,5,50);
En_Orig = 1./N.^(dim/2)*log2(N);


m = length(dim);
n = length(alpha);

En_Scaled = zeros(m,n);

for i = 1:m
    En_Scaled(i,:) = 1/N^(dim(i)/2)./alpha.^(dim(i)).*log2(alpha.^2*N);
end

% Variation of entropy under different scales for Original systems
figure;
plot(dim,En_Orig);
xlabel('Exponent');
ylabel('Entropy of original system');
grid on;
box on;




% Variation of entropy under different scales for scaled systems
figure;
surf(dim,alpha,En_Scaled,'EdgeColor','None');

xlabel('Exponent');
ylabel('Scales applied');
zlabel('Entropy of scaled system');

figure;
surf(dim,alpha,En_Scaled);
view(2);
xlabel('Exponent');
ylabel('Scales applied');