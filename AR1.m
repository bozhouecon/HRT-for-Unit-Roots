function Y = AR1(increments,mu,delta,rho)

T    = size(increments,1);
reps = size(increments,2);
Y    = zeros(T+1,reps);
for t=2:T+1
    Y(t,1:reps) = delta+rho*Y(t-1,1:reps)+increments(t-1,1:reps);
end
Y    = repmat(mu,T+1,reps) + Y;

end