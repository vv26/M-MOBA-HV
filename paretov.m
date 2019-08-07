function Vc=paretov(xb,sig,f1,r,n,tau,refpts)
xref=refpts(1);yref=refpts(2);xmin=refpts(3);ymin=refpts(4);
Vc=zeros(1,r);
for k=1:r
    nu=n(k)-1;
    tCDFc=gamma((nu+1)/2)/(sqrt(nu*pi)*gamma(nu/2));
    if f1(k)==1
        f2=zeros(1,r);
        for i=[1:k-1,k+1:r]
            a=xb(i,1)<=xb([1:k-1,k+1:r],1);
            b=xb(i,2)<=xb([1:k-1,k+1:r],2);
            c=a|b;
            cn=sum(c);
            f2(i)=(cn==(r-1));
        end
    else
        f2=f1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%here
    A=[xmin;sort(xb(f1==1,1));xref];
    B=[ymin;sort(xb(f1==1,2));yref];
    
    %%%%%%%%%%%%%%contribution of the zones on the left of f1 to the Vc of k-th point %%%%%%%%%%%%
    l=sum(f1);
    A2=A;
    B2=B;
    xband=[xb(k,1),min(A(A>xb(k,1)))];
    yband=[xb(k,2),min(B(B>xb(k,2)))];
    xyband=[xband,yband];
    for lj=1:l+1 %top lj-th layer
        yab=[B2(l-lj+2),B2(l-lj+3)];
        for li=1:lj %li is row
            xab=[A2(li),A2(li+1)];
            Vc(k)=Vc(k)+tdv(f1,f2,xab,yab,xb,sig,k,n,tau,1,tCDFc,refpts,xyband);           
        end
    end
    
    
    if f1(k)==1
        f4=(f2~=f1);
        f4(k)=0;
        %xband=[xb(k,1),min(A(A>xb(k,1)))];
        %yband=[xb(k,2),min(B(B>xb(k,2)))];
        %xyband=[xband,yband];
        A2=[xband(1);sort(xb(f4==1,1));xband(2)];
        B2=[yband(1);sort(xb(f4==1,2));yband(2)];
        
        %%%%%%%%%%%%%%contribution of the zones between f2 and f1 to the Vc of k-th point%%%%%%%%%%%%
        q=sum(f4);
        for lj=1:q+1 
            yab=[B2(q-lj+2),B2(q-lj+3)];
            for li=1:lj %li is row
                xab=[A2(li),A2(li+1)];
                Vc(k)=Vc(k)+tdv(f1,f2,xab,yab,xb,sig,k,n,tau,-1,tCDFc,refpts,xyband);
            end;
        end;
        
        %%%%%%%%%%%%%%contribution of the zones on the right of f2 to the Vc of k-th point %%%%%%%%%%%%
        p=sum(f2);
        A2=[sort(xb(f2==1,1));xref];
        B2=[sort(xb(f2==1,2));yref];
        pj=0;
        for li=1:p  
            yab=[B2(p-li+1),B2(p+1)];
            xab=[A2(li),A2(li+1)];
            pj=pj+td(xab,yab,xb,sig,k,n,tau);
        end;
        HVnc=hypervolumerr(xband,yband,f4,xb);
        Vc(k)=Vc(k)+pj*HVnc;   
    end
    
    
end;