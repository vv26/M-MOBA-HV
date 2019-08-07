function HVc4=hypervolumerr(xband,yband,f4,xb)%%%%%%to calculate the reduced HV when the k-th point moved to RHS of f2 pareto
HVc4=0;
q=sum(f4);
A=[xband(2);sort(xb(f4==1,1),'descend');xband(1)];
B=[yband(1);sort(xb(f4==1,2));yband(2)];
for i=1:q+1
    HVc4=HVc4+(B(i+1)-B(1))*(A(i)-A(i+1));
end
if (B(i+1)-B(1))*(A(i)-A(i+1))<0
    fprintf('error1');
end