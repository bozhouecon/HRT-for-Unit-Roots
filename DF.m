function   reject = DF(Y,statistic,model)

n = length(Y)-1;
z = Y(2:n+1,1);
switch statistic
    case 't-statistic'
        switch model
            case 'AR'
                X      = Y(1:n,1);
                b      = X\z;
                e      = z-X*b;
                covm   = ((e'*e)/n)/(X'*X);
                s      = sqrt(covm(1,1));
                DF     = (b(1,1)-1)/s;
                cv     = -1.95;
                reject = (DF<=cv);
            case 'ARD'
                X      = [Y(1:n,1) ones(n,1)];
                b      = X\z;
                e      = z-X*b;
                covm   = ((e'*e)/n)/(X'*X);
                s      = sqrt(covm(1,1));
                DF     = (b(1,1)-1)/s;
%                 cv     = -2.89; % n100
                cv     = -2.86; % n2500
                reject = (DF<=cv);
            case 'TS'
                X      = [Y(1:n,1) ones(n,1) (1:n)'];
                b      = X\z;
                e      = z-X*b;
                covm   = ((e'*e)/n)/(X'*X);
                s      = sqrt(covm(1,1));
                DF     = (b(1,1)-1)/s;
                cv     = -3.45;
                reject = (DF<=cv);
        end
    case 'estimator'
        switch model
            case 'AR'
                X      = Y(1:n,1);
                b      = X\z;
                DF     = n*(b-1);
                cv     = -7.8643;
                reject = (DF<=cv);
            case 'ARD'
                X      = [ones(n,1) Y(1:n,1)];
                b      = X\z;
                DF     = n*(b(2,1)-1);
%                 cv     = -13.52; % n100
                cv     = -14.05; % n2500
                reject = (DF<=cv);
            case 'TS'
                X      = [ones(n,1) Y(1:n,1) (1:n)'];
                b      = X\z;
                DF     = n*(b(2,1)-1);
                cv     = -20.4737;
                reject = (DF<=cv);
        end
end

end
