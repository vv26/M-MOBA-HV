function [xb0,sig0]=initialxs(r,n0,Mu,tsig)
xb0=zeros(r,2);sig0=zeros(r,2);
for i=1:n0
    X(1:r,2*i-1:2*i)= (normrnd(Mu,tsig));
    xb0= (xb0+X(1:r,2*i-1:2*i));
end
xb0 = (xb0/n0);
for i=1:n0
    sig0 = (sig0+(X(1:r,2*i-1:2*i)-xb0).^2);
end
sig0 = (sqrt(sig0/(n0-1)));