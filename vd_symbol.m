function f = vd_symbol(t,x)
        mu = 0.00001;
        dxdt_1 = [0 1];
        dxdt_2 = [-1 1-x(2)^2]./mu;
        f = [dxdt_1; dxdt_2];
                
end
