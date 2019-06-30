function  [Ig,score_R] = rankscores_approximate(R,g)

T = size(R,1);

if strcmp(g,'Gaussian')   
    score_R = norminv(R/(T+1),0,1);
    Ig      = 1; 
elseif strcmp(g,'DE')
    b_g     = 1/sqrt(2);
    score_R = sign(R/(T+1)-0.5)/b_g;
    Ig      = 2;
elseif strcmp(g,'t3')
    t_inv   = tinv(R/(T+1),3)/sqrt(3);
    score_R = 4*t_inv./(1+t_inv.^2);
    Ig      = 2;
elseif strcmp(g,'Logistic') 
    score_R = (pi/sqrt(3))*(1-(T+1-R)./R)./(1+(T+1-R)./R);
    Ig      = 1.0966;
else
    error('Distribution not yet implemented.')
end

end
