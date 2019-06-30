function [reject] = ERS_arma(Y,c_bar)

T             = length(Y);
alpha_bar     = 1+c_bar/T;
DeltaY        = Y(2:T)-Y(1:T-1);
S(1)          = sum(DeltaY'*DeltaY);
Ya(1,1)       = Y(1);
Ya(2:T,1)     = Y(2:T)-alpha_bar*Y(1:T-1);
Za(1,1)       = 1;
Za(2:T,1)     = 1-alpha_bar;
b             = (Za'*Za)\(Za'*Ya);
e             = Ya-Za*b;
S(2)          = sum(e'*e);

Dterm         = [DeltaY(8:end-1) DeltaY(7:end-2) DeltaY(6:end-3) DeltaY(5:end-4) DeltaY(4:end-5) DeltaY(3:end-6) DeltaY(2:end-7) DeltaY(1:end-8)];          
X             = [Y(9:end-1) Dterm ones(T-9,1)];
beta_hat      = (X'*X)\(X'*DeltaY(9:end));
e             = DeltaY(9:end)-X*beta_hat;
sigma2        = e'*e/T;   

omega2        = sigma2/((1-sum(beta_hat(2:end-1)))^2);
TS            = (S(2)-alpha_bar*S(1))/omega2;

if T<150
    cv = 3.11; 
elseif (150<=T)&&(T<=250)
    cv = 3.17; 
elseif T>250
    cv = 3.26;
else
    error('Choose another n...')
end

reject        = TS<cv;

end


