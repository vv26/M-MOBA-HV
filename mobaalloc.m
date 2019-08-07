function [recordmobaf1,recordmobaxb,recordmobasig,Vcr,ct,recordmobaequal,record2]=mobaalloc(xb0,sig0,n0,r,jn,budgets,budgeti,sps,Mu,f0,refpts,matlabzero)

xb=xb0; sig=sig0;n=n0*ones(r,1);budget=0;
nonalgorithm=0;
for k=1:r
    a=xb(k,1)<xb([1:k-1,k+1:r],1);
    b=xb(k,2)<xb([1:k-1,k+1:r],2);
    c=a|b;
    cn=sum(c);
    f1(k)=(cn==(r-1));%%if point k is on the pareto front based on five observations, f1(k)=1,0 otherwise
end;
for j=1:jn
    tbudget=budgets+(j-1)*budgeti;
    while(budget<tbudget)
        pj= (paretot(xb,sig,f1,r,n,1));
        
        if sum(pj)>matlabzero
            [m,mn]=max(pj);%%m is the largest number of pj and mn is the correspoding alternative's sequence number
            %[m,mn]= min(pj);
            recordmobaequal(j)=mn+100;
            %T(3,j)=T(3,j)+1;%%how many times use algorithm
        else
            pj=(paretot(xb,sig,f1,r,n,10));
            if sum(pj)>matlabzero
                [m,mn]=max(pj);
                recordmobaequal(j)=mn+200;
                %T(3,j)=T(3,j)+1;
            else
                mn=mod(nonalgorithm,r)+1;
                recordmobaequal(j)=mn+300;
                nonalgorithm=nonalgorithm+1;%% switch twice to tau=10
                %T(4,j)=T(4,j)+1;
            end
        end
        X=sps(mn,:,n(mn)-4);%++++++++++++++++++++++++++++++++++++++
        sig(mn,1:2)=sqrt((n(mn)-1)/n(mn)*sig(mn,1:2).^2+1/(n(mn)+1)*(X-xb(mn,1:2)).^2);
        xb(mn,1:2)=(n(mn)*xb(mn,1:2)+X)/(n(mn)+1);
        n(mn)=n(mn)+1;
        budget=budget+1;
        
        for k=1:r
            a=xb(k,1)<xb([1:k-1,k+1:r],1);
            b=xb(k,2)<xb([1:k-1,k+1:r],2);
            c=a|b;
            cn=sum(c);
            f1(k)=(cn==(r-1));
        end;
        record2(:,j)=pj;%%record in each iteration which alternative gets the sampling and the pj accordingly
    end
    
    recordmobaf1(:,j)=f1;
    recordmobaxb(:,1:2,j)=xb;
    recordmobasig(:,1:2,j)=sig;
    
    
    Vcr(j)=realVc(f1,f0,xb,Mu,refpts);
    ct(j)=(sum(f1==f0)==r);
end