function HVc=hypervolumer(f1,f2,xb,k,xab,yab,refpts,xyband)%%change of hypervolume when (x,y)between f2 and f1
HVc=zeros(1,4);
%xref=refpts(1);yref=refpts(2);xmin=refpts(3);ymin=refpts(4);
x1=sum(xab)/2;
y1=sum(yab)/2;
f4=(f2~=f1);
f4(k)=0;
%A=xb(f1==1,1);
%B=xb(f1==1,2);
%A=[xmin;sort(A);xref];
%B=[yref;sort(B,'descend');ymin];
xband=xyband(1:2);
yband=xyband(3:4);
A=[xband(1);sort(xb(f4==1,1));xband(2)];
B=[yband(2);sort(xb(f4==1,2),'descend');yband(1)];
f5=(A>x1)&(B>y1);
xm=min(A(B<y1));
ym=min(B(A<x1));
A=[A(f5);xm];
B=[ym;B(f5)];
q=sum(f5);
%HV=HV+(A(1)-x)*(B(1)-y);
% if ((A(1)-x1)*(B(1)-y1))<0
%     fprintf('errorR1\n');
% end
HVc(1)=1;HVc(2)=-B(1);HVc(3)=-A(1);HVc(4)=A(1)*B(1);
for i=2:q+1
    %HV=HV+(A(i)-A(i-1))*(B(i)-y);
    HVc(4)=HVc(4)+(A(i)-A(i-1))*B(i);
    HVc(3)=HVc(3)-(A(i)-A(i-1));
%     if ((A(i)-A(i-1))*(B(i)-y1))<0
%         fprintf('errorR2\n');
%     end
end
%HV=hypervolumerr(xband,yband,f4,xb)-HV;
HVc(1)=-HVc(1);HVc(2)=-HVc(2);HVc(3)=-HVc(3);HVc(4)=hypervolumerr(xband,yband,f4,xb)-HVc(4);
% if (HVc(1)*x1*y1+HVc(2)*x1+HVc(3)*y1+HVc(4))<-10^(-5)
%     fprintf('errorR3 ,%f\n',HVc(1)*x1*y1+HVc(2)*x1+HVc(3)*y1+HVc(4));
% end