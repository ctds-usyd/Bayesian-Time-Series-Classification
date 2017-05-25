
 % generated xvar and y data ------------------- 
clear all
close all
nind=28;
nobs=12;
mu_beta1_true=[5,10,0]';
mu_beta2_true=[5,-10,10]';
tau_beta1_true=[0.5,0.5,0.5]';
tau_beta2_true=[0.1,0.1,0.1]';
pi1_true=.5;
pi2_true=.5;
sigsq1_true=0.5;
sigsq2_true=0.5;
xvar=linspace(1./nobs,1,nobs)';
%xvar=linspace(1,nobs,nobs)';
for j=1:nind
    if(j<15)
        beta_true(:,j)=normrnd(mu_beta1_true,sqrt(tau_beta1_true));
        funct_true(:,j)=beta_true(1,j)+beta_true(2,j)*xvar+beta_true(3,j)*xvar.^2;
        y(:,j)=funct_true(:,j)+normrnd(zeros(nobs,1),sqrt(sigsq1_true)*ones(nobs,1));
        
    else
         beta_true(:,j)=normrnd(mu_beta2_true,sqrt(tau_beta2_true));
        funct_true(:,j)=beta_true(1,j)+beta_true(2,j)*xvar+beta_true(3,j)*xvar.^2;
        y(:,j)=funct_true(:,j)+normrnd(zeros(nobs,1),sqrt(sigsq2_true)*ones(nobs,1));
                
    end
end
figure
plot(xvar,y(:,1:14),'r')
hold
plot(xvar,y(:,16:28),'b')