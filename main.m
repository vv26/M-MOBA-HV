clear;%% use same random number and switch to equal allocation when Vj=0
n0=5;
matlabzero=10^(-6);%%%%the actual 0 in the program
% Mu=[5 1;4 2;1 5;6 0.5;7 0.25];
Mu=[1 5;5 1;3.2 2.1;3 2;2 3.1;6 4;5 5;4 6];
% Mu=[1 5;5 1;3 3;3.1 2;2 3.1;4 2.1;2.1 4;5.5 5;3.5 5;6 6];
%Mu=[1 5;2 4;4 2;5 1];
%Mu=[0.5 5.5;1.9 4.2;2.8 3.3;3 3;3.9 2.1;4.3 1.8;4.6 1.5;3.8 6.3;4.8 5.5;5.2 5;5.9 4.1;6.3 3.8;6.7 7.2;7 7;7.9 6.1;9 9];
% Mu=[1 8;2 5;3.5 5.01;3 2;2.5 8;3 7;3.05 2.2;1.5 6;1.05 8.03;2.1 5.2;2.5 4;2.6 3.9];
%Mu = [1 7;2 6.8;3 6;4 4;5 5];
[r,cl]=size(Mu);tsig=1*ones(r,2);
refpts=[10,10,-1000,-1000];%%%%%[xref yref xmin ymin]for hypervolume 
for k=1:r 
    a=Mu(k,1)<Mu([1:k-1,k+1:r],1);
    b=Mu(k,2)<Mu([1:k-1,k+1:r],2);
    c=a|b;
    cn=sum(c);
    f0(k)=(cn==(r-1));%% if point k is on the pareto front, f0(k)=1,0 otherwise
end

jn=400;   %%how many different values of budget will be tested
mre=2; %% the maximum repitation

ct=zeros(3,jn);%%the times of the right choice;ct(1,) equal allocation; ct(2,) allocation accoding to probability
Vcr=zeros(3,jn);%%the real change of hypervolume after j-th budget is allocated
%T=zeros(4,jn);%%how many times algorithm used with each budget.1, 2 for HV, 3, 4 for mmoba
budgets=0;%%the smallest budget
budgeti=1;%%budget step

record1=zeros(r,mre,jn);%% record in each iteration which alternative gets the sampling and the Vj accordingly
record2=zeros(r,mre,jn);%% record in each iteration which alternative gets the sampling and the pj accordingly
recordequalf1=zeros(r,mre,jn);%% record in equal allocation each iteration the selected pareto front
recordhvf1=zeros(r,mre,jn);%% record in hv allocation each iteration the selected pareto front
recordmobaf1=zeros(r,mre,jn);%% record in mmonba each iteration the selected pareto front
recordequalxb=zeros(r,2,mre,jn);
recordhvxb=zeros(r,2,mre,jn);
recordhvsig=zeros(r,2,mre,jn);
recordmobaxb=zeros(r,2,mre,jn);
recordmobasig=zeros(r,2,mre,jn);
recordmobaequal=zeros(mre,jn);
recordhvequal=zeros(mre,jn);

for re=1:mre
    sps=samplespace(Mu,tsig,jn,r,cl);
    %%%%%%%%%%%%%%%the initial sample mean value and variance%%%%%%%%%%%%%%%%%%
    [xb0,sig0]=initialxs(r,n0,Mu,tsig);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%euqal allocation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [recordequalf1(:,re,:),recordequalxb(:,:,re,:),Vcr1(re,:),ct1(re,:)]=equalalloc(xb0,n0,r,jn,budgets,budgeti,sps,Mu,f0,refpts);
% %%%%%%%%%%%%%%%% M-MOBA HV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [recordhvf1(:,re,:),recordhvxb(:,:,re,:),recordhvsig(:,:,re,:),Vcr2(re,:),ct2(re,:),recordhvequal(re,:),record1(:,re,:)]=hvalloc(xb0,sig0,n0,r,jn,budgets,budgeti,sps,Mu,f0,refpts,matlabzero);
    
% %%%%%%%%%%%%%%%%M-MOBA PCS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [recordmobaf1(:,re,:),recordmobaxb(:,:,re,:),recordmobasig(:,:,re,:),Vcr3(re,:),ct3(re,:),recordmobaequal(re,:),record2(:,re,:)]=mobaalloc(xb0,sig0,n0,r,jn,budgets,budgeti,sps,Mu,f0,refpts,matlabzero);
    fprintf('re=%d\n',re);
end
    
    
    
 for j=1:jn
    ct(1,j)=sum(ct1(:,j));
    ct(2,j)=sum(ct2(:,j));
    ct(3,j)=sum(ct3(:,j));

    Vcr(1,j)=sum(Vcr1(:,j));
    Vcr(2,j)=sum(Vcr2(:,j));
    Vcr(3,j)=sum(Vcr3(:,j));
 end
x=budgets+[0:(jn-1)]*budgeti;
% 
cp=ct/mre;
% 
Vcr_average=Vcr/mre;
% 
figure
plot(x,cp(1,1:jn),'r',x,cp(2,1:jn),'-g',x,cp(3,1:jn),'--b')
legend('Equal','M-MOBA HV','M-MOBA PCS')
xlabel('budget')
ylabel('P(CS)')


figure
plot(x,Vcr_average(1,1:jn),'r',x,Vcr_average(2,1:jn),'-g',x,Vcr_average(3,1:jn),'--b')
legend('Equal','M-MOBA HV','M-MOBA PCS')
xlabel('budget')
ylabel('change of hypervolume')
