function Vcr=realVc(f1,f0,xb,Mu,refpts)
xref=refpts(1);yref=refpts(2);
A=[Mu(f0==1,1);xb(f1==1,1)];
B=[Mu(f0==1,2);xb(f1==1,2)];
m0=sum(f0);m1=sum(f1);
m=m0+m1;

for k=1:m0
    a=A(k)<A([1:k-1,k+1:m]);
    b=B(k)<B([1:k-1,k+1:m]);
    c=a|b;
    cn=sum(c); 
    f5(k)=(cn==m-1);
end
q=sum(f5);

for k=(m0+1):m
    a=A(k)<A([1:k-1,k+1:m]);
    b=B(k)<B([1:k-1,k+1:m]);
    c=a|b;
    cn=sum(c); 
    f5(k)=(cn==m-1);
end;
p=sum(f5);

A1=A(f5);B1=B(f5);
Cf=[zeros(q,1);ones(p-q,1)];
bi=1:p;

for i=1:p-1
    for j=i+1:p
        if A1(bi(i))>A1(bi(j))
            temp=bi(i);
            bi(i)=bi(j);
            bi(j)=temp;
        end
    end
end

k=0;
for i=1:p-1
    if Cf(bi(i))~=Cf(bi(i+1))
        k=k+1;
        C1(k,1:2)=[A1(bi(i+1)),B1(bi(i))];
    end
end

m3=m-p+k;
        

A1=[sort(A1);xref];B1=[yref;sort(B1,'descend')];


V1=0;
for i=1:p
    V1=V1+(xref-A1(i))*(B1(i)-B1(i+1));
%     if (xref-A1(i))*(B1(i)-B1(i+1))<0
%         fprintf('errorVcr1\n');
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if k>0
    A2=sort([A(~f5);C1(:,1)]);B2=sort([B(~f5);C1(:,2)],'descend');
else
    A2=sort(A(~f5));B2=sort(B(~f5),'descend');
end
%%%%%%%%%%%%%%%%%test%%%%%%%%%%%%%%%%%%%%%%
[test1,test2]=size(A2);
if m3~=test1
     fprintf('m3_error\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A2=[A2;xref];B2=[yref;B2];
V2=0;
for i=1:m3
    V2=V2+(xref-A2(i))*(B2(i)-B2(i+1));
%     if (xref-A2(i))*(B2(i)-B2(i+1))<0
%         fprintf('errorVcr2\n');
%     end
end
Vcr=V1-V2;
if Vcr<0;
     fprintf('errorVcr3\n');
end



