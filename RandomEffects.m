% by Sally Cripps (2017)
nloop=250;
nwarmup=200;
ysort=y;
np=3;
ncomp=2;
%set hyperparameters and limits
nobs=12;
nind=28;
%xvar=linspace(1/nobs,1,nobs)';
%PUT IN SENSIBLE HYPERPARAMS
tau_prior_a=-1;
tau_prior_b=0;
tau_up_limit_beta=150;
%create storage space 
Beta_all=zeros(np,nind,nloop+nwarmup);
Sigsq_all=zeros(ncomp,nloop+nwarmup);
Mu_Beta_all=zeros(np,ncomp,nwarmup+nloop);
Tau_Beta_all=zeros(np,ncomp,nwarmup+nloop);
Y_fit_all=zeros(nobs,nind,nloop+nwarmup);
pi_comp1_all=zeros(nloop+nwarmup,1);
pcomp1_all=zeros(nind,nloop+nwarmup);
comp_all=zeros(nind,nloop+nwarmup);
%initialise parameters
comp=ones(nind,1);
 k1=linspace(1,2,2)';
 k2=linspace(3,28,26)';
u=rand(nind,1);
k1=find(u<0.5);
k2=find(u>0.5);
 comp(k2)=2;
n1=length(k1);
n2=nind-n1;
Beta=zeros(nind,np);
Tau_Beta=ones(np,ncomp);
Mu_Beta=ones(np,ncomp);
zmat=[ones(nobs,1) xvar xvar.^2];
zzmat=repmat(zmat,nind,1);
ystack=reshape(ysort,nobs*nind,1);
Sig_up_lim=var(ystack);
zdashz=zmat'*zmat;
zzinv=inv(zdashz);
ydev=zeros(nobs,nind);
for j=1:nind
    Beta(j,:)=zzinv*zmat'*ysort(:,j);
    ydev(:,j)=ysort(:,j)-zmat*Beta(j,:)';
end
tau_up_limit_beta=range(Beta).^2'/4;
mu_up_limit=max(Beta)';
mu_low_limit=min(Beta)';
Sigsq(1)=var(reshape(ydev(:,k1),nobs*n1,1));
Sigsq(2)=var(reshape(ydev(:,k2),nobs*n2,1));
Mu_Beta(:,1)=mean(Beta(k1,:));
Mu_Beta(:,2)=mean(Beta(k2,:));
Tau_Beta(:,1)=var(Beta(k1,:));
Tau_Beta(:,2)=var(Beta(k2,:));
pi_comp1=0.5;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     GIBBS LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
for p=1:nwarmup + nloop%For 1
   
    for j=1:nind
        if comp(j)==1
            prior_beta_mean=Mu_Beta(:,1);
            prior_beta_prec=diag(1./Tau_Beta(:,1));
            like_beta_prec=zdashz/Sigsq(1);
            like_beta_mean=zzinv*zmat'*ysort(:,j);
            post_beta_prec=like_beta_prec+prior_beta_prec;    
        else
            prior_beta_mean=Mu_Beta(:,2);
            prior_beta_prec=diag(1./Tau_Beta(:,2));
            like_beta_prec=zdashz/Sigsq(2);
            like_beta_mean=zzinv*zmat'*ysort(:,j);
            post_beta_prec=like_beta_prec+prior_beta_prec;
        end
        post_beta_var=post_beta_prec\eye(np);
        post_beta_var=0.5*(post_beta_var+post_beta_var');
        post_beta_mean=post_beta_var*(like_beta_prec*like_beta_mean+prior_beta_prec*prior_beta_mean);
        Beta(j,:)=mvnrnd(post_beta_mean,post_beta_var);
        ydev(:,j)=ysort(:,j)-zmat*Beta(j,:)';
    end
   
   if p>nwarmup
       Y_fit_all(:,:,p)=zmat*Beta';
   end
  % Drawing Sigmasq for group1
    k1=find(comp==1);
    n1=length(k1);
    ystack=reshape(ydev(:,k1),nobs*n1,1);
    sigsq_gampar_b=sum(ystack.^2)/2;
    sigsq_gampar_a=n1*nobs/2-1;
    Sigsq(1)=1/gamrnd(sigsq_gampar_a,1/sigsq_gampar_b);
     % Drawing Sigmasq for group2
    k2=find(comp==2);
    n2=length(k2);
    ystack=reshape(ydev(:,k2),nobs*n2,1);
    sigsq_gampar_b=sum(ystack.^2)/2;
    sigsq_gampar_a=n2*nobs/2-1;
    Sigsq(2)=1/gamrnd(sigsq_gampar_a,1/sigsq_gampar_b);
      
   %Drawing random effects parameters
   %First mean and variance parameters of intercepts
   %and positive slopes
   
   %Draw mean and variance of random effects model for first group
    if n1>0
        if n1==1
          Mu_Beta(:,1)=normrnd(Beta(k1,:)',sqrt(Tau_Beta(:,1)/n1)); 
        else
          Mu_Beta(:,1)=normrnd(mean(Beta(k1,:))',sqrt(Tau_Beta(:,1)/n1));
        end
    else
        Mu_Beta(:,1)=rand(np,1).*(mu_up_limit-mu_low_limit)+mu_low_limit; 
    end
   
    tau_beta_a=(n1/2-1)*ones(np,1);
    if n1>1
        tau_beta_b=sum((ones(n1,1)*Mu_Beta(:,1)'-Beta(k1,:)).^2)/2;
    else
        tau_beta_b=(ones(n1,1)*Mu_Beta(:,1)'-Beta(k1,:)).^2/2;
    end
    if n1>2
        u=rand*ones(np,1);
        const1=gamcdf(1./tau_up_limit_beta,tau_beta_a,1./tau_beta_b');
        const2=ones(np,1)-u.*(1-const1);
        Tau_Beta(:,1)=1./gaminv(const2,tau_beta_a,1./tau_beta_b');
    else
       tau_temp=rand(np,1).*tau_up_limit_beta;
       if n1==0
           met_rat=1;
       else
           tau_loglike_new=sum(-(1+tau_beta_a).*log(tau_temp)-tau_beta_b'./tau_temp);
           tau_loglike_old=sum(-(1+tau_beta_a).*log(Tau_Beta(:,1))-tau_beta_b'./Tau_Beta(:,1));
           met_rat=min(1,exp(tau_loglike_new-tau_loglike_old));
       end
       u=rand;
       if u<met_rat
         Tau_Beta(:,1)=tau_temp;
       end
    end
   %Draw mean and variance of random effects model for second  group
   if n2>0
        if n2==1
          Mu_Beta(:,2)=normrnd(Beta(k2,:)',sqrt(Tau_Beta(:,2)/n2)); 
        else
          Mu_Beta(:,2)=normrnd(mean(Beta(k2,:))',sqrt(Tau_Beta(:,2)/n2));
        end
    else
        Mu_Beta(:,2)=rand(np,1).*(mu_up_limit-mu_low_limit)+mu_low_limit; 
    end
   
    tau_beta_a=(n2/2-1)*ones(np,1);
    if n2>1
        tau_beta_b=sum((ones(n2,1)*Mu_Beta(:,2)'-Beta(k2,:)).^2)/2;
    else
        tau_beta_b=(ones(n2,1)*Mu_Beta(:,2)'-Beta(k2,:)).^2/2;
    end
    if n2>2
        u=rand*ones(np,1);
        const1=gamcdf(1./tau_up_limit_beta,tau_beta_a,1./tau_beta_b');
        const2=ones(np,1)-u.*(1-const1);
        Tau_Beta(:,2)=1./gaminv(const2,tau_beta_a,1./tau_beta_b');
    else
       tau_temp=rand(np,1).*tau_up_limit_beta;
       if n2==0
           met_rat=1;
       else
           tau_loglike_new=sum(-(1+tau_beta_a).*log(tau_temp)-tau_beta_b'./tau_temp);
           tau_loglike_old=sum(-(1+tau_beta_a).*log(Tau_Beta(:,2))-tau_beta_b'./Tau_Beta(:,2));
           met_rat=min(1,exp(tau_loglike_new-tau_loglike_old));
       end
       u=rand;
       if u<met_rat
         Tau_Beta(:,2)=tau_temp;
       end
    end
   
   %Allocating indidivuals to  Components
        for j=1:nind
            %Likelihood of compoenent 1 after integrating beta out
            A1=zmat'*ysort(:,j)/Sigsq(1)+diag(1./Tau_Beta(:,1))*Mu_Beta(:,1);
            V1=eye(np)/(zdashz/Sigsq(1)+diag(1./Tau_Beta(:,1)));
            quadform=A1'*V1*A1;
            sumysq=sum(ysort(:,j).^2)/Sigsq(1);
            summusq=Mu_Beta(:,1)'*diag(1./Tau_Beta(:,1))*Mu_Beta(:,1);
            like_comp1=(2*pi)^(-nobs/2)*Sigsq(1)^(-nobs/2)*prod(Tau_Beta(:,1))^-0.5*exp(-.5*(sumysq+summusq-quadform))*det(V1);
             %Likelihood of component 2 after integrating beta out
             A2=zmat'*ysort(:,j)/Sigsq(2)+diag(1./Tau_Beta(:,2))*Mu_Beta(:,2);
            V2=eye(np)/(zdashz/Sigsq(2)+diag(1./Tau_Beta(:,2)));
            quadform=A2'*V2*A2;
            sumysq=sum(ysort(:,j).^2)/Sigsq(2);
            summusq=Mu_Beta(:,2)'*diag(1./Tau_Beta(:,2))*Mu_Beta(:,2);
            like_comp2=(2*pi)^(-nobs/2)*Sigsq(2)^(-nobs/2)*prod(Tau_Beta(:,2))^-0.5*exp(-.5*(sumysq+summusq-quadform))*det(V2);
            pcomp1(j)=like_comp1*pi_comp1/(like_comp1*pi_comp1+like_comp2*(1-pi_comp1));
            u=rand;
            if u<pcomp1(j)
                comp(j)=1;
                else
                comp(j)=2;
            end
        end
        %Drawing prob of components
        k1=find(comp==1);
        n1=length(k1);
        pi_comp1=betarnd(n1+1,nind-n1+1);
        pi_comp2=1-pi_comp1;
      
        
     %Storing Parameter Values
    Beta_all(:,:,p)=Beta';
    Sigsq_all(:,p)=Sigsq;
    Mu_Beta_all(:,:,p)=Mu_Beta;
    Tau_Beta_all(:,:,p)=Tau_Beta;
    pi_comp1_all(p)=pi_comp1;
    pcomp1_all(:,p)=pcomp1;
    comp_all(:,p)=comp;
    
end
 
%End for 1
%end 
% see each individual class
comp
for j=1:n1
    figure
hold
title('Fit')
  plot(xvar,mean(Y_fit_all(:,j,nwarmup+1:nloop),3),'b')  
  plot(xvar,ysort(:,j),'*')  
end


for j=n1+1:28
   figure
hold
title('Fit ')  
  plot(xvar,mean(Y_fit_all(:,j,nwarmup+1:nloop),3),'r')  
 plot(xvar,ysort(:,j),'*')  
end