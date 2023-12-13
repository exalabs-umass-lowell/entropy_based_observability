% Reduce 10 by 10 matrix A to 5 by 5 with different MOR methods
% Full order of system is not observable, while reduced order system is
% observable.

% Brandon.yang---2016

%   A    square matrix (10 by 10)
%   B    vector     (1 by 10)
%   C    vecotr     (1 by 10)



A = [-1 0 3 0 0 0 0 0 0 -0.1695;
    -0.0410 0 3 0 0 0 0 0 0 0;
    0 -314 0 0 0 0 0 0 0 0;
    0 0 -0.866 -20 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    -0.9386 -51.9849 -0.7999 0 0 -41.1153 0 0 0 0;
    0 0 0 0 0 0 0 -5 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0];


B = zeros(1,10);
B(1,2) = 0.0926;
B(1,6) = 0.4428;
B(1,7) = 2.1179;
B(1,8) = 2.1179;
B = B';
C = zeros(1,10);
C(1,1) = 0.4777;
C(1,3) = -0.3239;


% Based on the rank of observability matrix to check if it's observable
Ob = obsv(A,C);
Of = rank(Ob);        % rank of the full order system
Q = Ob*Ob;
Metric = min(eig(Q))/max(eig(Q));




% Arnoldi iteration method
[Q,H] = Arnoldi(A,B,4);          % (A,B,N) -- N is the reduced order by Arnoldi iteration method

AR = Q'*A*Q;       % order_reduced A matrix
BR = Q'*B;         % order_reduced B matrix
CR = C*Q;          % order_reduced C matrix
sys1 = ss(A,B,C,0);     % original system
sys_reduced = ss(AR,BR,CR,0);   % order_reduced system


Obr = obsv(AR,CR);
Or = rank(Obr);    % rank of the reduced order system
Qr = Obr*Obr;
Metricr = min(eig(Qr))/max(eig(Qr));



