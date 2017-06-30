function y = rs(t,x)

    %parameters for the model
    a = 0.2; 
    b = 0.2; 
    c = 5.7;       
%%% equations

y = [
     -x(2)-x(3);            % dx/dt
     x(1)+a*x(2);           % dy/dt
     b+x(3)*(x(1)-c);       % dz/dt
] ; 