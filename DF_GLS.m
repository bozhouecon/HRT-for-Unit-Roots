function reject = DF_GLS(Y,model,c_bar)

T         = length(Y);
alpha_bar = 1+c_bar/T;
switch model
    case 'ARD'
        Z  = ones(T,1); 
        cv = -1.95;
    case 'TS'
        Z  = [ones(T,1) (1:1:T)'];      
end
y_a      = [Y(1);Y(2:T)-alpha_bar*Y(1:T-1)];
Z_a      = [Z(1,:);Z(2:T,:)-alpha_bar*Z(1:T-1,:)];
beta_hat = Z_a\y_a;
yd       = Y - Z*beta_hat;

X      = yd(1:T-1,1);
z      = yd(2:T,1);
b2     = X\z;
e      = z-X*b2;
covm   = ((e'*e)/(T-1))/(X'*X);
s      = sqrt(covm(1,1));
DF     = (b2-1)/s;

reject = DF<=cv;

end