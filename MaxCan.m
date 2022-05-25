function [Val,M] = MaxCan(Cannonic,Y)

[M,I] = max(Cannonic);

switch I
    
    case 1
        
        Val = Y(1,1);
        
    case 2
        
        Val = Y(1,2);
        
    case 3
        
        Val = Y(1,3);
        
    case 4
        
        Val = Y(1,4);
        
    case 5
        
        Val = Y(1,5);
        
    case 6
        
        Val = Y(1,6);
    
end

end