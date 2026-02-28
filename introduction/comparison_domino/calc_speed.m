function v = calc_speed(t,phi,n,d)
% Calculation of intrinsic collision speed
c = 2; % average over c dominos
threshold = 1e-3;

v = zeros(n,1);
v(1:3) = ones(3,1)*NaN;
for i = 4:n
    % first index where |phi_i| exceeds threshold
    ind1 = find(abs(phi(i, :)) > threshold, 1, 'first');

    % first index where |phi_(i-c)| exceeds threshold
    ind2 = find(abs(phi(i-c, :)) > threshold, 1, 'first');
    Delta_t = t(ind1) - t(ind2);
    Delta_x = d*c;
    v(i) = Delta_x/Delta_t;
end
plot(1:n,v,'o-')
xlabel('domino sequence number')
ylabel('intrinsic collision speed')