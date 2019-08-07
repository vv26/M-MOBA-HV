function HVc=hypervolumel(f1,f2,xb,k,xab,yab,refpts,xyband)%%change of hypervolume when (x,y)on the left of f1
xref=refpts(1);yref=refpts(2);xmin=refpts(3);ymin=refpts(4);
x1=sum(xab)/2;
y1=sum(yab)/2;
A=xb(f1==1,1);
B=xb(f1==1,2);
A=[xmin;sort(A);xref];
B=[yref;sort(B,'descend');ymin];
fv=(A>x1)&(B>y1);
VXm=min(A(B<y1));
VYm=min(B(A<x1));
VX=[A(fv);VXm];
VY=[VYm;B(fv)];
p=sum(fv)+1;

HVc(1)=1;HVc(2)=-VY(1);HVc(3)=-VX(1);HVc(4)=VX(1)*VY(1);
%HV=(VX(1)-x)*(VY(1)-y);
for i=2:p
    HVc(3)=HVc(3)-(VX(i)-VX(i-1));
    HVc(4)=HVc(4)+(VX(i)-VX(i-1))*(VY(i));
    %HV=HV+(VX(i)-VX(i-1))*(VY(i)-y);
end


if f1(k)==1
        xband=xyband(1:2);
        yband=xyband(3:4);
        if (y1<yband(2))&&(y1>yband(1))
            f4=(f2~=f1);
            f4(k)=0;
            q=sum(f4);
            A=[xband(1);sort(xb(f4==1,1));xband(2)];
            B=[yband(2);sort(xb(f4==1,2),'descend');yband(1)];
            for i=1:q+1
                if y1>B(i+1)
                    HVc(3)=HVc(3)+A(i+1)-A(1);
                    HVc(4)=HVc(4)-(A(i+1)-A(1))*B(i+1);
                    %HV=HV+(A(i+1)-A(1))*(y-B(i+1));
                    for j=i+1:q+1
                        HVc(4)=HVc(4)+(A(j+1)-A(1))*(B(j)-B(j+1));
                        %HV=HV+(A(j+1)-A(1))*(B(j)-B(j+1));
                    end
                    break;
                end
            end
        elseif (x1<xband(2))&&(x1>xband(1))
            f4=(f2~=f1);
            f4(k)=0;
            q=sum(f4);
            A=[xband(2);sort(xb(f4==1,1),'descend');xband(1)];
            B=[yband(1);sort(xb(f4==1,2));yband(2)];
            for i=1:q+1
                if x1>A(i+1)
                    HVc(2)=HVc(2)+B(i+1)-B(1);
                    HVc(4)=HVc(4)-(B(i+1)-B(1))*A(i+1);
                    %HV=HV+(B(i+1)-B(1))*(x-A(i+1));
                    for j=i+1:q+1
                        HVc(4)=HVc(4)+(B(j+1)-B(1))*(A(j)-A(j+1));
                        %HV=HV+(B(j+1)-B(1))*(A(j)-A(j+1));
                    end
                    break;
                end
            end
        elseif (y1>yband(2))||(x1>xband(2))
            f4=(f2~=f1);
            f4(k)=0;
            HVc(4)=HVc(4)+hypervolumerr(xband,yband,f4,xb);
        end    
end