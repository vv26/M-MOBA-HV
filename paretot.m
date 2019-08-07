function pj=paretot(xb,sig,f1,r,n,tau)
for k=1:r
    pj(k)=0;
    if f1(k)==1
        f2=zeros(1,r);
        for i=[1:k-1,k+1:r]
            a=xb(i,1)<=xb([1:k-1,k+1:r],1);
            b=xb(i,2)<=xb([1:k-1,k+1:r],2);
            c=a|b;
            cn=sum(c);
            f2(i)=(cn==(r-1));
        end;
    else
        f2=f1;
    end;
    A=sort(xb(f2==1,1));
    B=sort(xb(f2==1,2));
    l=sum(f2);

    if f1(k)==1
        if l==sum(f1)-1
            A1=[-inf;A;inf];
            B1=[-inf;B;inf];
            for i=1:l+1
                xab=[A1(i),A1(i+1)];yab=[B1(l-i+2),B1(l-i+3)];
                pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau);
            end;
        else
            A2=[-inf;A;inf];
            B2=[-inf;B;inf];
            xab(1)=max(A2(A2<xb(k,1)));yab(1)=max(B2(B2<xb(k,2)));
            xab(2)=min(A2(A2>xb(k,1)));yab(2)=min(B2(B2>xb(k,2)));
            pj(k)=td(xab,yab,xb,sig,k,n,tau);
        end;
    else
        A=[A;inf];
        for i=1:l
            xab=[A(i),A(i+1)];yab=[B(l-i+1),inf];
            pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau);
        end;
    end; 
end;
pj=1-pj;