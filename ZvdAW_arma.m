function [reject] = ZvdAW_arma(Y,g,c_bar,AR)

    Delta_Y     = diff(Y); 
    
    if  strcmp(AR,'AR1') 
        theta          = regress(Delta_Y(2:end),Delta_Y(1:end-1)); 
        eps_hat        = Delta_Y(2:end)-theta*Delta_Y(1:end-1); 
    elseif strcmp(AR,'AR2')
        theta          = regress(Delta_Y(3:end),[Delta_Y(2:end-1) Delta_Y(1:end-2)]); 
        eps_hat        = Delta_Y(3:end)-theta(1,1)*Delta_Y(2:end-1)-theta(2,1)*Delta_Y(1:end-2); 
    elseif strcmp(AR,'AR3')
        theta          = regress(Delta_Y(4:end),[Delta_Y(3:end-1) Delta_Y(2:end-2) Delta_Y(1:end-3)]); 
        eps_hat        = Delta_Y(4:end)-theta(1,1)*Delta_Y(3:end-1)-theta(2,1)*Delta_Y(2:end-2)-theta(3,1)*Delta_Y(1:end-3);   
    elseif strcmp(AR,'AR4')
        theta          = regress(Delta_Y(5:end),[Delta_Y(4:end-1) Delta_Y(3:end-2) Delta_Y(2:end-3) Delta_Y(1:end-4)]); 
        eps_hat        = Delta_Y(5:end)-theta(1,1)*Delta_Y(4:end-1)-theta(2,1)*Delta_Y(3:end-2)-theta(3,1)*Delta_Y(2:end-3)-theta(4,1)*Delta_Y(1:end-4);   
    elseif strcmp(AR,'AR8')
        theta          = regress(Delta_Y(9:end),[Delta_Y(8:end-1) Delta_Y(7:end-2) Delta_Y(6:end-3) Delta_Y(5:end-4) Delta_Y(4:end-5) Delta_Y(3:end-6) Delta_Y(2:end-7) Delta_Y(1:end-8)]); 
        eps_hat        = Delta_Y(9:end)-theta(1,1)*Delta_Y(8:end-1)-theta(2,1)*Delta_Y(7:end-2)-theta(3,1)*Delta_Y(6:end-3)-theta(4,1)*Delta_Y(5:end-4)-theta(5,1)*Delta_Y(4:end-5)-theta(6,1)*Delta_Y(3:end-6)-theta(7,1)*Delta_Y(2:end-7)-theta(8,1)*Delta_Y(1:end-8);  
    else
        error('Not an available choice yet.')
    end
    
    R           = tiedrank(eps_hat);
    sigma_hat   = sqrt(var(eps_hat(:,1)));
    T           = length(R);
    
    if  strcmp(g,'Gaussian') 
        sigma_eg_hat   = cov(eps_hat(:,1),norminv(R(:,1)/(T+1),0,1));
        sigma_eg_hat   = sigma_eg_hat(1,2);
        rho            = sigma_eg_hat/sigma_hat;
        cv             = 0.9610+1.8836*rho-3.9772*rho^2+6.7373*rho^3-5.4544*rho^4+1.6916*rho^5;
    elseif strcmp(g,'DE') 
        sigma_eg_hat   = cov(eps_hat(:,1),(sign(R(:,1)/(T+1)-0.5)*sqrt(2)));
        sigma_eg_hat   = sigma_eg_hat(1,2);
        rho            = sigma_eg_hat/sigma_hat;
        cv             = 0.2467+2.3014*rho-3.5841*rho^2+4.2983*rho^3-2.4461*rho^4+0.5398*rho^5;
    elseif strcmp(g,'t3') 
        t_inv          = tinv(R(:,1)/(T+1),3)/sqrt(3);
        sigma_eg_hat   = cov(eps_hat(:,1),4*t_inv./(1+t_inv.^2));
        sigma_eg_hat   = sigma_eg_hat(1,2);
        rho            = sigma_eg_hat/sigma_hat;
        cv             = 0.2467+2.3014*rho-3.5841*rho^2+4.2983*rho^3-2.4461*rho^4+0.5398*rho^5;
    else
        error('Distribution not yet implemented.')
    end
    
    LR          = LLR(eps_hat,R,c_bar,g,sigma_hat,rho);
    reject      = LR>cv;
    
end


