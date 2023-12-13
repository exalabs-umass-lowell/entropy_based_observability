x0 = 0;
x1 = -2;


% Establish the van der pol model
t  = 0:0.05:100;   % time scale
xa = [x0 x1];
[t,x] = ode45('vd', t, xa); 

plot(x(:,1),x(:,2));
xlabel('Simulation Time(s)');
ylabel('State');
zlabel('observability index');
%title('');
legend('State 1', 'State 2');
grid on;

% m = x(101:400,:);
% a = unique(m(:,1));
% out1 = [a,histc(x(:),a)];
% M = sum(out1(:,2));
% Freq1 = out1(:,2)/M;
% plot(out1(:,1),Freq1);


% 
% 
% u=linspace(-15, 15, 100);
% q=quantizer([9 6], 'float');
% y=quantize(q, u);
% plot(u, y); title(tostring(q))


%%
t = [4.95:.05:20]; % Times at which to sample the sine function
sig = x(99:400,1); % Original signal, a sine wave
partition = [-1.9:0.1:2]; % 300 intervals
codebook = [-2.0:.1:2]; % 301

[index,quants] = quantiz(sig,partition,codebook); % Quantize.
% plot(t,sig,'x',t,quants,'.')
% grid on;
% xlabel('Simulation Time(s)');
% ylabel('State Variation');

%hold on;

sig2 = x(99:400,2); % Original signal, a sine wave
[index,quants2] = quantiz(sig2,partition,codebook); % Quantize.
% plot(t,sig2,'o',t,quants2,'--');
% grid on;
% xlabel('Simulation Time(s)');
% ylabel('State Variation');
% legend('Original signal-2','Quantized signal-2');





hold on;
plot(quants,quants2);
hold on;
plot(x(101:400,1),x(101:400,2));
xlabel('State x1');
ylabel('State x2');
legend('Quantized signal','Original signal');

%%




%%








% 
% 
% 
% 
% b = unique(m(:,2));
% out2 = [b,histc(x(:),b)];
% 
% 
% % M = size(x,1);
% % % Probability distribution
% % 
% % p = 
% % Entro = sum(p.*log(1./p));
% 
% 
% 
% 
% 
% a = unique(x);
% out = [a,histc(x(:),a)];