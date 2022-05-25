function FB = FilterBank(S,Nf,LP,HP,Fn)

[m,n] = size(S);

for i = 1:m
    
    for j = 1:length(LP)
        
        FB(i,1:n,j) = BandPass(S(i,:),Nf,LP(1,j),HP(1,j),Fn);
        
    end
    
end

end