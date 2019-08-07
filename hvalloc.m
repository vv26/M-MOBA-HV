function [recordhvf1,recordhvxb,recordhvsig,Vcr,ct,recordhvequal,record1]=hvalloc(xb0,sig0,n0,r,jn,budgets,budgeti,sps,Mu,f0,refpts,matlabzero)

xb=xb0; sig=sig0;n=n0*ones(r,1);budget=0;%T=zeros(2,jn);
nonalgorithm=0;%% to record how many times will we have to use equal allocation

for k=1:r
    a=xb(k,1)<xb([1:k-1,k+1:r],1);
    b=xb(k,2)<xb([1:k-1,k+1:r],2);
    c=a|b;
    cn=sum(c);
    f1(k)=(cn==(r-1));%%if point k is on the pareto front based on five observations, f1(k)=1,0 otherwise
end
for j=1:jn
    %fprintf('re=%d ',re);
    %fprintf('eoc j=%d\n',j);
    tbudget=budgets+(j-1)*budgeti;
    while(budget<tbudget)
        Vc= paretov(xb,sig,f1,r,n,1,refpts);
        %%%%%%%%%%%%%%%%to make sure Vc>0%%%%%%%%%%%%%%%%
%         if sum(Vc<(-matlabzero))>0
%             fprintf('error VC<0,when j=%d, ',j);
%             for testi=1:r
%                 %                     if Vc(testi)<0
%                 fprintf('%f ,',Vc(testi));
%                 %                     end
%             end
%             fprintf('\n');
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if abs(sum(Vc))>matlabzero
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [m,mn]=max(Vc);%%m is the largest number of Vj and mn is the correspoding alternative's sequence number
            %[m,mn]= min(Vj);
            recordhvequal(j)=mn+100;
            %T(1,j)=T(1,j)+1;%%how many times use algorithm
        else
            Vc=paretov(xb,sig,f1,r,n,10,refpts);
            %%%%%%%%%%%%%%%%to make sure Vc>0%%%%%%%%%%%%%%%%
%             if sum(Vc<(-matlabzero))>0
%                 fprintf('error VC<0,when j=%d, ',j);
%                 for testi=1:r
%                     %                         if Vc(testi)<0
%                     fprintf('%f ,',Vc(testi));
%                     %                         end
%                 end
%                 fprintf('\n');
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if abs(sum(Vc))>matlabzero
                [m,mn]=max(Vc);
                recordhvequal(j)=mn+200;
                %T(1,j)=T(1,j)+1;
            else
                %T(2,j)=T(2,j)+1;
                mn=mod(nonalgorithm,r)+1;
                recordhvequal(j)=mn+300;
                nonalgorithm=nonalgorithm+1;
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
        end
        record1(:,j)=Vc;%%record in each iteration which alternative gets the sampling and the Vj accordingly
    end
    recordhvf1(:,j)=f1;
    recordhvxb(:,1:2,j)=xb;
    recordhvsig(:,1:2,j)=sig;
    
    
    Vcr(j)=realVc(f1,f0,xb,Mu,refpts);
    ct(j)=(sum(f1==f0)==r);
    
end