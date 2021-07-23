 function dxdt=rhs(t,x)
        dxdt = x/(0.01-x^2);