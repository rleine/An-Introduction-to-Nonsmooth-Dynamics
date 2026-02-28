function prox = prox_CT(z,a)
% proximal point on C_{Ti} = [-a_i, a_i], elementwise
    prox = max(-a,min(z,a));
end