function y = Cubic_Root(a,c,d)

q = -d/2/a;
r = c/3/a;
y = (q+sqrt(q^2+r^3))^(1/3) + (q+sqrt(q^2+r^3))^(1/3);

end