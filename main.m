%% 
clear 
clc
%% Generate increments
warmup  = 50;
T       = 100;   
reps    = 20000;  
f       = 't3';        % Gaussian t1 t2 t3 DE Logistic Pearson stable
sigma   = 1;
theta   = -0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: IID
epsilons = innovation(T+warmup,reps,f);
increments = sigma*epsilons(warmup+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: Heteroskedasticity
% % Model = garch('Constant',1,'GARCH',0.65,'ARCH',0.25,'Distribution',struct('Name','t','DoF',3));
% Model = garch('Constant',1,'GARCH',0.19,'ARCH',0.8);
% for i = 1:reps
%     [~,epsilons(:,i)] = simulate(Model,T+warmup);
% end
% increments = sigma*epsilons(warmup+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: SerialCorreclation MA
% v(1,:)      = sigma*epsilons(1,:);
% for i = 2:T+warmup
%    v(i,:)      = sigma*epsilons(i,:)+theta*sigma*epsilons(i-1,:);
% end
% increments = v(warmup+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: SerialCorreclation AR
% v(1,:)      = sigma*epsilons(1,:);
% for i = 2:T+warmup
%    v(i,:)      = theta*v(i-1,:)+sigma*epsilons(i,:);
% end
% increments = v(warmup+1:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DGP: SerialCorreclation ARMA
% v(1,:)      = sigma*epsilons(1,:);
% for i = 2:T+warmup
%    v(i,:)      = theta*v(i-1,:) + sigma*epsilons(i,:)+theta*sigma*epsilons(i-1,:);
% end
% increments = v(warmup+1:end,:);

%% Test Power
c = -0:-0.5:-30;
g = 'DE';  % Gaussian DE t3
for i = 1:length(c)
    i
    beta    = 1+c(i)/T;
%     beta    = 1+c(i)/(n-n*tau);
    Y       = AR1(increments,1,0,beta); % AR1(innovations,mu,delta,rho) 
for s=1:reps
    [ZvdAW_rej(s)]     = ZvdAW(Y(2:end,s),g,-7);
%     [ZvdAW_rej(s)]     = ZvdAW_symmetric(Y(2:end,s),g,-7);
%     [ZvdAW_rej(s)]     = ZvdAW_sc_homos(Y(2:end,s),g,-7);
% %     [ZvdAW_rej(s)]     = ZvdAW_fhat_split(Y(2:end,s),-7,tau,0.35);
% %     [Jansson_rej(s)]   = Jansson_fhat_split(Y(2:end,s),-7,tau,0.35); 
%     [ZvdAW_rej(s)]     = ZvdAW_fhat(Y(2:end,s),-7,0.35);
    [Jansson_rej(s)]   = Jansson_fhat(Y(2:end,s),-7,0.5); 
    [ERS_rej(s)]       = ERS(Y(2:end,s),-7);
%     [DF_GLS_rej(s)]    = DF_GLS(Y(2:end,s),'ARD',-7);
%     [DF_rho_rej(s)]    = DF(Y(2:end,s),'estimator','ARD');
%     [DF_t_rej(s)]      = DF(Y(2:end,s),'t-statistic','ARD');
end
    POWER_ZvdAW(i)   = mean(ZvdAW_rej);
    POWER_Jansson(i) = mean(Jansson_rej);
    POWER_ERS(i)     = mean(ERS_rej);
%     POWER_DF_rho(i)  = mean(DF_rho_rej);
%     POWER_DF_t(i)    = mean(DF_t_rej);
%     POWER_DF_GLS(i)  = mean(DF_GLS_rej);
end
%% Plot
load('PowerEnvelope_Gaussian.mat');
load('PowerEnvelope_DE_semi.mat');
load('PowerEnvelope_Jf2.mat')
figure;
% plot(-c,PowerEnvelope_Jf2_semi(1:length(c)),'k-','LineWidth',1.5);hold on;
% plot(-c,PowerEnvelope_Gaussian(1:length(c)),'k-','LineWidth',1.5);hold on;
plot(-c,POWER_ZvdAW(1:length(c)),'r:');hold on;
plot(-c,POWER_Jansson(1:length(c)),'b.--');hold on;
plot(-c,POWER_ERS(1:length(c)),'g--');hold on;
% plot(-c,POWER_DF_rho(1:length(c)),'r.--');hold on;
% plot(-c,POWER_DF_t(1:length(c)),'r--');hold on;
% plot(-c,POWER_DF_GLS(1:length(c)),'r:');hold on;
% legend('POWERENVELOPE-Gaussian','POWERENVELOPE-DE-semi','ZvdAW','DF-rho','DF-t','DF-GLS','ERS','Location','southeast')
% text(0,POWER_ZvdAW(1,1),['size:' num2str(POWER_ZvdAW(1,1))])



