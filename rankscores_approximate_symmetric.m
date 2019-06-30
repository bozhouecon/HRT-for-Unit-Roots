function  [Ig,epsilon_R] = rankscores_approximate_symmetric(R_plus,signs,g)

T = size(R_plus,1);
    
if strcmp(g,'Gaussian')   
    epsilon_R = signs.*norminv((T+1+R_plus)/(2*(T+1)),0,1);
    Ig        = 1; 
elseif strcmp(g,'DE')
    epsilon_R = signs.*sign((T+1+R_plus)/(2*(T+1))-0.5)*sqrt(2);
    Ig        = 2;
elseif strcmp(g,'t3')
    t_inv     = tinv((T+1+R_plus)/(2*(T+1)),3)/sqrt(3);
    epsilon_R = signs.*4.*t_inv./(1+t_inv.^2);
    Ig        = 2;
else
    error('Distribution not yet implemented.')
end

end
