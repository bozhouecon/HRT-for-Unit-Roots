function [LR,score_R] = LLR_symmetric(epsilon_hat,R,signs,c_bar,g,sigma_hat)

T            = size(epsilon_hat,1);
[Ig,score_R] = rankscores_approximate_symmetric(R,signs,g);
Zg           = cumsum([zeros(1,1);score_R])/sqrt(T);
Z            = [zeros(1,1);cumsum(epsilon_hat)/sqrt(T)]/sigma_hat;  
LR           = sum(Z.*[diff(Zg);zeros(1,1)])*c_bar-0.5*Ig*mean(Z.^2)*c_bar^2;

end