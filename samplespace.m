function sps=samplespace(Mu,tsig,jn,r,cl)
for i=1:r
    for j=1:cl
        sps(i,j,:)=normrnd(Mu(i,j),tsig(i,j),jn,1);
    end;
end;