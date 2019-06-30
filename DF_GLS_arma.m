function reject = DF_GLS_arma(Y,model,c_bar)

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
Yd       = Y - Z*beta_hat;

y      = Yd(10:end);
DeltaY = diff(Yd); 
Dterm  = [DeltaY(8:end-1) DeltaY(7:end-2) DeltaY(6:end-3) DeltaY(5:end-4) DeltaY(4:end-5) DeltaY(3:end-6) DeltaY(2:end-7) DeltaY(1:end-8)];          
X      = [Yd(9:end-1) Dterm];
b      = X\y;
e      = y-X*b;
covm   = ((e'*e)/(T-1))/(X'*X);
s      = sqrt(covm(1,1));
DF     = (b(1)-1)/s;

reject = DF<=cv;

end