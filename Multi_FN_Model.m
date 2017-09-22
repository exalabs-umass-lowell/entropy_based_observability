% This function is about four connected F
function y = Multi_FN_Model(t,x)

v1 = x(1);
r1 = x(2);
v2 = x(3);
r2 = x(4);
v3 = x(5);
r3 = x(6);
v4 = x(7);
r4 = x(8);
c = -1/0.08;
a = -0.7;
b = -0.8;
I = 0; % Step input
k = 0.4; % Coupling strength



y = [(v1-v1^3/3+r1+I+k*(v2-v1));
    -1/c*(v1-a+b*r1)+10;
    (v2-v2^3/3+r2+k*(v1-v2)+k*(v3-v2));
    -1/2/c*(v2-a+b*r2);
    (v3-v3^3/3+r3+k*(v2-v3)+k*(v4-v3));
    -1/5/c*(v3-a+b*r3);
    (v4-v4^3/3+r4+k*(v3-v4));
    -1/10/c*(v4-a+b*r4)];
end