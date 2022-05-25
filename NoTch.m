function BP = NoTch(S,Nf,LP,HP,Fn)

[z,p,k] = butter(Nf,[LP HP]/(Fn/2),'stop');
sos = zp2sos(z,p,k);

%fvt = fvtool(sos,'Fs',Fn);

VF(:,1) = S;

for i = 1:Nf
    
    VF(:,i + 1) = filter(sos(i,1:3),sos(i,4:6),VF(:,i));
    
end

BP = VF(:,Nf + 1);

end