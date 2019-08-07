
CH=zeros(1,r);
CM=zeros(1,r);
for i=1:mre
    for j=1:jn
        n=mod(recordhvequal(i,j),100);
        if n~=0%%%%%when j=1, tbudget=0, no alternative gets a new sample. at this time, n=0.
            CH(n)=CH(n)+1;
        end    
        n=mod(recordmobaequal(i,j),100);
        if n~=0
            CM(n)=CM(n)+1;
        end
        
    end
end
CM=CM/(jn*mre);
CH=CH/(jn*mre);

        