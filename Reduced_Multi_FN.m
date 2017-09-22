

function y = Reduced_Multi_FN(t,x)

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

N = x(7);

if N == 7
    y = [-c*(v1-v1^3/3+r1+I+k*(v2-v1));
        (v1-a+b*r1)-c*10;
        -2*c*(v2-v2^3/3+r2+k*(v1-v2)+k*(v3-v2));
        (v2-a+b*r2);
        -5*c*(v3-v3^3/3+r3+k*(v2-v3)+k*(Cubic_Root(10*c/3,-10*c*(1-k),-10*c*(r4))-v3));
        (v3-a+b*r3);
        0;
        (Cubic_Root(10*c/3,-10*c*(1-k),-10*c*(r4))-a+b*r4)];
    
else if N==6
        y = [-c*(v1-v1^3/3+r1+I+k*(v2-v1));
            (v1-a+b*r1)-c*10;
            -2*c*(v2-v2^3/3+r2+k*(v1-v2)+k*(Cubic_Root(10*c/3,-10*c*(1-k)-v2,-10*c*(r3))));
            (v2-a+b*r2);
            0;
            (Cubic_Root(5*c/3,-5*c*(1-k),-5*c*(r3))-a+b*r3);
            0;
            (Cubic_Root(10*c/3,-10*c*(1-k),-10*c*(r4))-a+b*r4)];
    else if N == 5
            y = [-c*(v1-v1^3/3+r1+I+k*(Cubic_Root(10*c/3,-10*c*(1-k),-10*c*(r2))-v1));
                (v1-a+b*r1)-c*10;
                0;
                (Cubic_Root(2*c/3,-2*c*(1-k),-2*c*(r2))-a+b*r2);
                0;
                (Cubic_Root(5*c/3,-5*c*(1-k),-5*c*(r3))-a+b*r3);
                0;
                (Cubic_Root(10*c/3,-10*c*(1-k),-10*c*(r4))-a+b*r4)];
        else if N==4
                y = [0;
                    (Cubic_Root(c/3,-c*(1-k),-c*(r4)-a+b*r1))-c*10;
                    0;
                    (Cubic_Root(2*c/3,-2*c*(1-k),-2*c*(r2)-a+b*r2));
                    0;
                    (Cubic_Root(5*c/3,-5*c*(1-k),-5*c*(r3))-a+b*r3);
                    0;
                    (Cubic_Root(10*c/3,-10*c*(1-k),-10*c*(r4))-a+b*r4)];
            end
        end
    end
end
end