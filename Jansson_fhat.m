function [reject] = Jansson_fhat(Y,c_bar,h)

T         = length(Y);
Delta_Y   = diff(Y); 

beta      = Y(1:end-1)\Y(2:end);
eps_hat_f = Y(2:end)-beta*Y(1:end-1);
eps_hat   = Delta_Y(1:end);
sigma_hat = sqrt(var(eps_hat));

u          = (eps_hat*ones(1,T-1)-ones(T-1,1)*eps_hat_f')/h;
score_fhat = sum(exp(-0.5*u.^2).*u/h,2)./sum(exp(-0.5*u.^2),2);

Ig           = mean(score_fhat.^2);
Zg           = cumsum([0;score_fhat])/sqrt(T);
Z            = [zeros(1,1);cumsum(eps_hat)/sqrt(T)]/sigma_hat;  

LR           = (sum((Z-mean(Z)).*[diff(Zg);0])+Z(end,:).*mean(Z))*c_bar-0.5*(Ig*mean(Z.^2)-(Ig-1)*(mean(Z)).^2)*c_bar^2;
cv           = 1.5510+0.9716*Ig-0.8488*Ig^2+0.1767*Ig^3-0.0054*Ig^4-0.0022*Ig^5;

reject       = LR>cv;

end


