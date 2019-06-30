function [innovations]=innovation(n,reps,f)

if strcmp(f,'Gaussian')
    innovations = normrnd(0,1,n,reps);
elseif strcmp(f,'skewnormal')
    delta       = 0.975;
    skewness    = (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^(1.5);
    kurtosis    = 3+2*(pi-3)*(delta*sqrt(2/pi))^4/(1-2*delta^2/pi)^2;
    innovations = snrnd(0,1,delta,n,reps);
elseif strcmp(f,'skewt4')
    lambda    = -0.5;
    innovations = skewtrnd(4,lambda,n,reps);
elseif strcmp(f,'DE')
    b_f         = 1/sqrt(2);
    U           = rand(n,reps)-0.5;    
    innovations = -b_f*sign(U).*log(1-2*abs(U));
elseif strcmp(f,'Logistic')
    U=rand(n,reps);
    innovations = (sqrt(3)/pi)*log(U./(1-U));
elseif strcmp(f,'Pearson')
    innovations = pearsrnd(0,1,3,36,n,reps);
elseif strcmp(f,'stable')
    innovations = stblrnd(0.5,0,1,0,n,reps);
elseif strcmp(f,'NIG')
    innovations = nigrnd(0.5,0,0,1,n,reps);
elseif strcmp(f,'t1')
    innovations = trnd(1,n,reps);
elseif strcmp(f,'t2')
    innovations = trnd(2,n,reps);
elseif strcmp(f,'t3')
    innovations = trnd(3,n,reps)/sqrt(3);
elseif strcmp(f,'t4')
    innovations = trnd(4,n,reps);
elseif strcmp(f,'t5')
    innovations = trnd(5,n,reps);
elseif strcmp(f,'t6')
    innovations = trnd(6,n,reps);
elseif strcmp(f,'t15')
    innovations = trnd(15,n,reps);
elseif strcmp(f,'t60')
    innovations = trnd(60,n,reps);
elseif strcmp(f,'t250')
    innovations = trnd(250,n,reps);
else
    error('Distribution not yet implemented.')
end


