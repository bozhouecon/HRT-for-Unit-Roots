function [reject] = ZvdAW_symmetric(Y,g,c_bar)

    Delta_Y(:,1)        = Y(2:end) - Y(1:end-1);
    epsilon_hat         = Delta_Y;
    sigma_hat           = sqrt(var(epsilon_hat));
    R                   = tiedrank(abs(epsilon_hat));
    signs               = sign(epsilon_hat);
    [LR,score_R]        = LLR_symmetric(epsilon_hat,R,signs,c_bar,g,sigma_hat);
    if  strcmp(g,'Gaussian') 
        sigma_eg_hat        = cov(epsilon_hat,score_R);
        sigma_eg_hat        = sigma_eg_hat(1,2);
        rho                 = sigma_eg_hat/sigma_hat;
        cv                  = 0.2121+2.6525*rho-0.8845*rho^2-0.4740*rho^3+0.6595*rho^4-0.3249*rho^5;
    elseif strcmp(g,'DE') 
        sigma_eg_hat        = cov(epsilon_hat,score_R);
        sigma_eg_hat        = sigma_eg_hat(1,2);
        rho                 = sigma_eg_hat/sigma_hat;
        cv                  = -1.1246+2.9303*rho-0.8734*rho^2+0.4786*rho^3-0.3415*rho^4+0.0883*rho^5;
    elseif strcmp(g,'t3') 
        sigma_eg_hat        = cov(epsilon_hat,score_R);
        sigma_eg_hat        = sigma_eg_hat(1,2);
        rho                 = sigma_eg_hat/sigma_hat;
        cv                  = -1.1246+2.9303*rho-0.8734*rho^2+0.4786*rho^3-0.3415*rho^4+0.0883*rho^5;
    else
        error('Distribution not yet implemented.')
    end
    
    reject             = LR>cv;
    
end


