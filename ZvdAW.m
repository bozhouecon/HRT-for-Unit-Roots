function [reject] = ZvdAW(Y,g,c_bar)

    T           = length(Y);
    increments  = Y(2:end) - Y(1:end-1);
    eps_hat     = increments;
    sigma_hat   = sqrt(var(eps_hat));
    R           = tiedrank(eps_hat);  
    
    if  strcmp(g,'Gaussian') 
        sigma_eg_hat        = cov(eps_hat,norminv(R/T,0,1));
        sigma_eg_hat        = sigma_eg_hat(1,2);
        rho                 = sigma_eg_hat/sigma_hat;
        cv                  = 0.9610+1.8836*rho-3.9772*rho^2+6.7373*rho^3-5.4544*rho^4+1.6916*rho^5;
    elseif strcmp(g,'DE') 
        sigma_eg_hat        = cov(eps_hat,(sign(R/T-0.5)*sqrt(2)));
        sigma_eg_hat        = sigma_eg_hat(1,2);
        rho                 = sigma_eg_hat/sigma_hat;
        cv                  = 0.2467+2.3014*rho-3.5841*rho^2+4.2983*rho^3-2.4461*rho^4+0.5398*rho^5;
    elseif strcmp(g,'t3') 
        t_inv               = tinv(R/T,3)/sqrt(3);
        sigma_eg_hat        = cov(eps_hat,4*t_inv./(1+t_inv.^2));
        sigma_eg_hat        = sigma_eg_hat(1,2);
        rho                 = sigma_eg_hat/sigma_hat;
        cv                  = 0.2467+2.3014*rho-3.5841*rho^2+4.2983*rho^3-2.4461*rho^4+0.5398*rho^5;
    else
        error('Distribution not yet implemented.')
    end
    
    LR          = LLR(eps_hat,R,c_bar,g,sigma_hat,rho);
    reject      = LR>cv;
    
end


