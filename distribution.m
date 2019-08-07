H=zeros(r,mre);
M=zeros(r,mre);
for i=1:mre
    for j=1:jn
        n=mod(recordhvequal(i,j),100);
        if n~=0%%%%%when j=1, tbudget=0, no alternative gets a new sample. at this time, n=0.
            H(n,i)=H(n,i)+1;
        end 
        n=mod(recordmobaequal(i,j),100);
        if n~=0
            M(n,i)=M(n,i)+1;
        end
    end
end


figure

bar(H(1,:))

