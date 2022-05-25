function BP = BandPass(S,Nf,LP,HP,Fn)

[z,p,k] = butter(Nf,[LP HP]/(Fn/2),"Bandpass");
sos = zp2sos(z,p,k);

VF(:,1) = S;

for i = 1:Nf
    
    VF(:,i + 1) = filter(sos(i,1:3),sos(i,4:6),VF(:,i));
    
end

BP = VF(:,Nf + 1);

end