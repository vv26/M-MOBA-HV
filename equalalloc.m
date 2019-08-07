function [recordequalf1,recordequalxb,Vcr,ct]=equalalloc(xb0,n0,r,jn,budgets,budgeti,sps,Mu,f0,refpts)
xb=xb0;n=n0*ones(r,1);budget=0;

for j=1:jn
    tbudget=budgets+(j-1)*budgeti;
    while(budget<tbudget)
        mn=mod(budget,r)+1;
        budget=budget+1;
        X=sps(mn,:,n(mn)-4);%++++++++++++++++++4 need to change when n0!=5
        xb(mn,1:2)= ((n(mn)*xb(mn,1:2)+X)/(n(mn)+1));
        n(mn)=n(mn)+1;
    end;
    for k=1:r
        a=xb(k,1)<xb([1:k-1,k+1:r],1);
        b=xb(k,2)<xb([1:k-1,k+1:r],2);
        c=a|b;
        cn=sum(c);
        f1(k)=(cn==(r-1));
    end;
    recordequalf1(:,j)=f1;
    recordequalxb(:,1:2,j)=xb;
    Vcr(j)=realVc(f1,f0,xb,Mu,refpts);
    ct(j)=(sum(f1==f0)==r);
end;