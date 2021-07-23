function vel = get_vel(density, params)

qc = params.qc;
kj = params.kj;
vf = params.vf;
kc = params.kc;

w = qc/(kj-kc);

k = density;

if((k>=0) && (k <= kc))
    vel = vf;
end

if (k > kc) && (k <= kj)
    vel = (qc - w*(k-kc))/k;
end