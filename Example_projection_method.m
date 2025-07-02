A = [2 1;
     1 2];

b = [2;-2];

[x,y] = LCP_projection_method(A,b,b,.1,1e-8,1000)