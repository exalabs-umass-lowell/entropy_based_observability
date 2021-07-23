function k = get_congested_density(q,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
qc = params.qc;
kj = params.kj;
vf = params.vf;
kc = params.kc;

w = qc/(kj-kc);

k = kc + (1/w)*(qc-q);

end

