function Ws = Windowing(S,Wt,Tt)

Wf = Wt*(1/16);

Tf = Tt*(1/16);

contwh = 1;

for i = 0:length(S)*Tf:length(S)-length(S)*Wf
    
    Whz = zeros(1,length(S));
    Who = ones(1,(length(S)*Wf));
    
    Whz(1,i + 1:(length(S)*Wf) + i) = Who;
    
    Ws(contwh,1:length(S)) = S.*Whz';
    P = mean(Ws(contwh,i + 1:(length(S)*Wf) + i));
    Ws(contwh,i + 1:(length(S)*Wf) + i) = (Ws(contwh,i + 1:(length(S)*Wf) + i) - P)/max(abs(Ws(contwh,i + 1:(length(S)*Wf) + i)));
    
    contwh = contwh + 1;
    
end

end