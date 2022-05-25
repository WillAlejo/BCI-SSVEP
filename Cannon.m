function [RoWhS,RoWhC] = Cannon(X,f,t)

[m,n] = size(X);

for i = 1:m
    
    for j = 1:length(f)
        
        fx = f(1,j);
        
        SENO = sin(2*pi*fx*t);
        
        Y = SENO;
        
        [A,B,r] = canoncorr(X(i,:)', Y');
        
        RoS(1,j) = r;
        
    end
    
    for k = 1:length(f)
        
        fx = f(1,k);
        
        COSENO = cos(2*pi*fx*t);
        
        Y = COSENO;
        
        [A,B,r] = canoncorr(X(i,:)', Y');
        
        RoC(1,k) = r;
        
    end
    
    RoWhS(i,1:length(RoS)) = RoS;
    
    RoWhC(i,1:length(RoC)) = RoC;
    
end

end


