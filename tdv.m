function Vc=tdv(f1,f2,xab,yab,xb,sig,k,n,tau,lr,tCDFc,refpts,xyband)%%%HyperVolume change
kp=n(k)*(n(k)+tau)./(tau*sig(k,1:2).^2);
mu=xb(k,1:2);
nu=n(k)-1;

if lr==1
    HVc=hypervolumel(f1,f2,xb,k,xab,yab,refpts,xyband);%%hypervolume coefficients: Hvc(1)xy+Hvc(2)xy+Hvc(3)y+Hvc(4)
else
    HVc=hypervolumer(f1,f2,xb,k,xab,yab,refpts,xyband);
end

%%%%%%%%%%%%%%%%%%%%%%%%int of (x-mu1)*StPDF(without the constant coeffient StCdfc) from xab1 to xab2 (or (y-mu2)*StPDF from yab1 to yab2)%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vint1=nu/(kp(1)*(1-nu))*(1+kp(1)*(xab(2)-mu(1))^2/nu)^((-nu+1)/2)-nu/(kp(1)*(1-nu))*(1+kp(1)*(xab(1)-mu(1))^2/nu)^((-nu+1)/2);
Vint2=nu/(kp(2)*(1-nu))*(1+kp(2)*(yab(2)-mu(2))^2/nu)^((-nu+1)/2)-nu/(kp(2)*(1-nu))*(1+kp(2)*(yab(1)-mu(2))^2/nu)^((-nu+1)/2);

%%%%% Hvc(1)xy+Hvc(2)xy+Hvc(3)y+Hvc(4)=H1(x-mu1)(y-mu2)+H2(x-mu1)+H3(y-mu2)+H4
H1=HVc(1);
H2=HVc(2)+HVc(1)*mu(2);
H3=HVc(3)+HVc(1)*mu(1);
H4=HVc(1)*mu(1)*mu(2)+HVc(2)*mu(1)+HVc(3)*mu(2)+HVc(4);
StCdfc=sqrt(kp(1:2))*tCDFc;%%constant in cdf of St
txab = (xab-mu(1))*sqrt(kp(1));
tyab = (yab-mu(2))*sqrt(kp(2));
pj1 = tcdf(txab(2),nu)-tcdf(txab(1),nu);
pj2 = tcdf(tyab(2),nu)-tcdf(tyab(1),nu);
Vc=H1*StCdfc(1)*StCdfc(2)*Vint1*Vint2+H2*StCdfc(1)*Vint1*pj2+H3*StCdfc(2)*Vint2*pj1+H4*pj1*pj2;