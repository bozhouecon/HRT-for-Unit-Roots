function [reject] = ZvdAW_fhat_split(Y,c_bar,tau,h)

T         = length(Y);
Delta_Y   = diff(Y); 

n1        = T*tau;
beta      = Y(end-n1-1:end-1)\Y(end-n1:end);
eps_hat_f = Y(end-n1:end)-beta*Y(end-n1-1:end-1);
pts       = linspace(.01,.99,99);
Fy        = ksdensity(eps_hat_f, pts, 'function','icdf', 'width',h);
Finv      = @(u) interp1(pts,Fy,u,'linear','extrap');

n2        = T - n1 - 1;
eps_hat   = Delta_Y(1:end-n1);
sigma_hat = sqrt(var(eps_hat));
R         = tiedrank(eps_hat);  
eps_R     = Finv(R/(n2+1))*sigma_hat;

u            = (eps_R*ones(1,n1+1)-ones(n2,1)*eps_hat_f')/h;
score_R_fhat = sum(exp(-0.5*u.^2).*u/h,2)./sum(exp(-0.5*u.^2),2);
Bg           = cumsum([0;score_R_fhat])/sqrt(n2);
Z            = [zeros(1,1);cumsum(eps_hat)/sqrt(n2)]/sigma_hat; 
Ig           = mean((score_R_fhat).^2);

LR           = (sum(Z.*[diff(Bg);0])+Z(end,:).*mean(Z))*c_bar-0.5*(Ig*mean(Z.^2)-(Ig-1)*(mean(Z)).^2)*c_bar^2;
cv           = 1.5510+0.9716*Ig-0.8488*Ig^2+0.1767*Ig^3-0.0054*Ig^4-0.0022*Ig^5;

reject       = LR>cv;

end


