function y=Reduced_Neuron_Three_HR(t,X)

%Initial conditiosn and parameter values

a = 1;
b = 3.0;
c = 1.0;
d = 5.0;
s = 4.0;
r = 0.006;
x0 = -1.6;
Epsilon = 0.5;
%B3 = 1000;
B2 = 50;
B1 = 10;

Excit = 4;
I_ext = 3;

if size(X,1) == 8    %After SP for the 1st one
    x1 = X(1);
    y1 = X(2);
    z1 = X(3);
    
    x2 = X(4);
    y2 = X(5);
    z2 = X(6);
    
    x3 = X(7);
    z3 = X(8);   % states
    dx2 = B2*(y2-a*x2^3+b*x2^2-z2+I_ext+Epsilon*((x1-x2)+(x3-x2)));
    y = [B1*(y1-a*x1^3+b*x1^2-z1+I_ext+Epsilon*(x2-x1));
        B1*(c-d*x1^2-y1);
        (s*(x1-x0)-z1)+Excit;
        B2*(y2-a*x2^3+b*x2^2-z2+I_ext+Epsilon*((x1-x2)+(x3-x2)));
        B2*(c-d*x2^2-y2);
        (s*(x2-x0)-z2);
        (s*(x3-x0)-z3-Epsilon*dx2)/(2*(b-d)*x3-3*a*x3^2-Epsilon);
        (s*(x3-x0)-z3)];
    
else if size(X,1) == 7   %After SP for the 2nd neuron
        x1 = X(1);
        y1 = X(2);
        z1 = X(3);
        
        x2 = X(4);
        z2 = X(5);
        
        x3 = X(6);
        z3 = X(7);
        dx2 = B2*(c-d*x2^2-a*x2^3+b*x2^2-z2+I_ext+Epsilon*((x1-x2)+(x3-x2)));
        dx1 = B1*(y1-a*x1^3+b*x1^2-z1+I_ext+Epsilon*(x2-x1));
        y = [B1*(y1-a*x1^3+b*x1^2-z1+I_ext+Epsilon*(x2-x1));
            B1*(c-d*x1^2-y1);
            (s*(x1-x0)-z1)+Excit;
            (s*(x2-x0)-z2-Epsilon*dx1)/(2*(b-d)*x2-3*a*x2^2-2*Epsilon);
            (s*(x2-x0)-z2);
            (s*(x3-x0)-z3-Epsilon*dx2)/(2*(b-d)*x3-3*a*x3^2-Epsilon);
            (s*(x3-x0)-z3)];
        
    else if size(X,1) == 6  %After SP for the 3rd neuron
            x1 = X(1);
            z1 = X(3);
            
            x2 = X(2);
            z2 = X(4);
            
            x3 = X(5);
            z3 = X(6);
            
            deno = 1/(2*(b-d)*x1-3*a*x1^2-Epsilon);
            dx2 = B2*(c-d*x2^2-a*x2^3+b*x2^2-z2+I_ext+Epsilon*((x1-x2)+(x3-x2)));
            dx1 = B1*(c-d*x1^2-a*x1^3+b*x1^2-z1+I_ext+Epsilon*(x2-x1));            
            y = [(s*(x1-x0)-z1-Epsilon*dx2)*deno; 
                (s*(x2-x0)-z2-Epsilon*dx1)/(2*(b-d)*x2-3*a*x2^2-2*Epsilon);
                (s*(x1-x0)-z1)+Excit;
                (s*(x2-x0)-z2);
                (s*(x3-x0)-z3-Epsilon*dx2)/(2*(b-d)*x3-3*a*x3^2-Epsilon);
                (s*(x3-x0)-z3)];
            
        end
    end
    
    
end



end
