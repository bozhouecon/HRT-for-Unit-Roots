%% 
clear 
clc
%% Generate increments
warmup  = 50;
n       = 250;   
reps    = 20000;  
f       = 't3';        % Gaussian t1 t2 t3 DE Logistic Pearson stable
sigma   = 3;
theta   = -0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: IID
epsilons = innovation(n+warmup,reps,f);
% increments = sigma*epsilons(warmup+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: SerialCorreclation MA
% v(1,:)      = sigma*epsilons(1,:);
% for i = 2:n+warmup
%    v(i,:)      = sigma*epsilons(i,:)+theta*sigma*epsilons(i-1,:);
% end
% increments = v(warmup+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: SerialCorreclation AR
% v(1,:)      = sigma*epsilons(1,:);
% for i = 2:n+warmup
%    v(i,:)      = theta*v(i-1,:)+sigma*epsilons(i,:);
% end
% increments = v(warmup+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: SerialCorreclation ARMA
v(1,:)      = sigma*epsilons(1,:);
for i = 2:n+warmup
   v(i,:)      = theta*v(i-1,:) + sigma*epsilons(i,:)+theta*sigma*epsilons(i-1,:);
end
increments = v(warmup+1:end,:);

%% Test Power
c = 0:-0.5:-30; 
g = 'Gaussian';  % Gaussian DE t3
for i = 1:length(c)
    beta    = 1+c(i)/n;
    Y       = AR1(increments,1,0,beta); % AR1(innovations,mu,delta,rho) 
for s=1:reps
    [ZvdAW_rej(s)]     = ZvdAW_arma(Y(2:end,s),g,-7,'AR8');
    [ERS_rej(s)]       = ERS_arma(Y(2:end,s),-7);
%     [DF_GLS_rej(s)]    = DF_GLS_arma(Y(2:end,s),'ARD',-7);
%     [ADF_rho_rej(s)]   = adftest(Y(2:end,s),'model','ARD','test','t2','lags',8);
%     [ZvdAW_rej(s)]     = ZvdAW(Y(2:end,s),g,-7);
%     [DF_t_rej(s)]      = DF(Y(2:end,s),'t-statistic','ARD');
%     [DF_rej(s)]        = adftest(Y(2:end,s),'model','ARD');
end
%     POWER_ADF_rho(i) = mean(ADF_rho_rej);
    POWER_ERS(i)     = mean(ERS_rej);
%     POWER_DF_GLS(i)  = mean(DF_GLS_rej);
    POWER_ZvdAW(i)   = mean(ZvdAW_rej);
end
%% Plot 
load('PowerEnvelope_Gaussian.mat');
load('PowerEnvelope_Jf2.mat')
load('PowerEnvelope_Jf2_semi.mat')
figure;
% plot(-c,PowerEnvelope_Gaussian(1:length(c)),'k-','LineWidth',1.5);hold on;
% plot(-c,PowerEnvelope_Jf2_semi(1:length(c)),'k-','LineWidth',1.5);hold on;
% plot(-c,POWER_ADF_rho(1:length(c)),'r.--');hold on;
% plot(-c,POWER_DF_GLS(1:length(c)),'r:');hold on;
plot(-c,POWER_ERS(1:length(c)),'g.--');hold on;
plot(-c,POWER_ZvdAW(1:length(c)),'b.--');hold on;
% legend('POWERENVELOPE-Gaussian','POWERENVELOPE-DE-semi','ZvdAW','DF-rho','DF-t','DF-GLS','ERS','Location','southeast')
% text(0,POWER_ZvdAW(1,1),['size:' num2str(POWER_ZvdAW(1,1))])


