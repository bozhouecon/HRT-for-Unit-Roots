function LR = LLR(epsilon_hat,R,c_bar,g,sigma_hat,rho)

T            = size(epsilon_hat,1);
[Ig,score_R] = rankscores_approximate(R,g);
Bg           = cumsum([zeros(1,1);score_R])/sqrt(T);
Z            = [zeros(1,1);cumsum(epsilon_hat)/sqrt(T)./sigma_hat];  
LR           = (sum(Z.*[diff(Bg);zeros(1,1)])+Z(end,:).*mean(Z)*rho)*c_bar...
                -0.5*(Ig*mean(Z.^2)-(Ig-rho^2)*(mean(Z)).^2)*c_bar^2;
                  
end