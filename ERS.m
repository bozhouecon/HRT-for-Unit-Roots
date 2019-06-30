function [reject] = ERS(Y,c_bar)

T             = length(Y);
alpha_bar     = 1+c_bar/T;
DeltaX        = Y(2:T,1)-Y(1:T-1,1);
S(1,1)        = sum(DeltaX'*DeltaX);
Ya(1,1)       = Y(1);
Ya(2:T,1)     = (Y(2:T)-alpha_bar*Y(1:T-1));
Za(1,1)       = 1;
Za(2:T,1)     = 1-alpha_bar;
b             = (Ya'*Za)/(Za'*Za);
e             = Ya-Za*b;
S(2,1)        = sum(e'*e);

DeltaY        = Y(2:T,1)-Y(1:T-1,1);
Y             = [ones(T-1,1) Y(1:T-1)];
beta          = (Y'*Y)\(Y'*DeltaY);
e             = DeltaY-Y*beta;
hat_sigma2    = e'*e/T;   
TS            = (S(2)-alpha_bar*S(1))/(hat_sigma2);

if T<=150
    cv = 3.11; 
elseif (150<=T)&&(T<=250)
    cv = 3.17; 
elseif T>=250
    cv = 3.26;
else
    error('Choose another n...')
end

reject        = TS<cv;

end


