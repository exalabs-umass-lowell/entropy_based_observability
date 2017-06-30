function f = vd(t,x)
        mu = 0.001;
        dxdt_1 = x(2);
        dxdt_2 = (x(2) - x(1) - x(2)^3/3)/mu;
        f = [dxdt_1; dxdt_2];
                
end
