%% Characteristics
clear;
clc;

tic
for sub = 1:5
    
    %sub = 1;
    for ses = 1:2
        
        %ses = 1;
        %Load Database (specific values)
        load("Sub" + sub + "_" + ses + "_multitarget.mat")
        
        fn = Data.AmpSamlingFrequency;
        Sub = strcat("Sub",num2str(sub));
        for NS = 1:10
            
            %NS = 1;
            EEG = Data.EEG(:,NS);
            
            %Time and frecuency vectors
            t = 0 : 1/fn : ((length(EEG)) - 1)/fn;
            f = 0 : fn/length(EEG) : fn - (fn/length(EEG));
            
            %FFT
            Sfft = ((abs(fft(EEG))).^2)./(length(EEG));
            
            %Band Pass Filter from database (Function)
            EEGF = BandPass(EEG,6,5,50,fn);
            
            %FFT filtered database
            SfftF = ((abs(fft(EEGF))).^2)./(length(EEG));
            
            %Windowing signal (Function) Wds = [NWindow,EEGWindow]
            Wds = Windowing(EEGF,1,0.25);
            
            %Frecuencies from CannonCorr
            fx = [6 6.5 7 7.5 8.2 9.3];
            
            %Frecuencies from CannonCorr (armonics)
            for i = 1:6
                fxA(i,:) = [fx(i) fx(i)*2 fx(i)*3 fx(i)*4 fx(i)*5];
            end
            
            %Filter Bank (Function) EEGFB = [NWindow,EEGWindow,NFilterBank]
            LF = [4 10 16 22 28 34];
            HF = [48.5 48.5 48.5 48.5 48.5 48.5];
            
            EEGFB = FilterBank(Wds,6,LF,HF,fn);

            %CCR
            %{
            %for j = 1:6
                for i = 1:6
                    [CCSFBA(:,i),CCCFBA(:,i)] = Cannon(Wds(:,:),fx(1,i),t);
                end
            %end
            %}
            %CCR (Armonics)
            %{
            for j = 1:6
                %for i = 1:6
                    [CCSFBA(:,:,j),CCCFBA(:,:,j)] = Cannon(Wds(:,:),fxA(j,:),t);
                %end
            end
            %}
            %FBCCR (Armonics) Cannon[Sen,cos] =
            %{Sen[NWindow,NArmonic,NFilterBank,NNomFrec],Cos[NWindow,NArmonic,NFilterBank,NNomFrec]}
            %
            for j = 1:6
                for i = 1:6
                    [CCSFBA(:,:,i,j),CCCFBA(:,:,i,j)] = Cannon(EEGFB(:,:,i),fxA(j,:),t);
                end
            end
            %
            %Real mean CCR 
            %{
            [m1,n1] = size(CCSFBA);
            %for k = 1:p1
                %for j = 1:o1
                    %{
                    for i = 1:m1
                        CCSFBAW(i,j,k) = mean(CCSFBA(i,:,j,k));
                        CCCFBAW(i,j,k) = mean(CCCFBA(i,:,j,k));
                    end
                    %}
                    for i = 1:n1
                        NwRe = 1;
                        for Nw = 1:8:m1 - 4    %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
                            CCSFBAP(NwRe,i) = mean(CCSFBA(Nw:Nw + 4,i));   %Nw:Nw + 2,i,j,k      Nw:Nw + 4,i,j,k       Nw:Nw + 4,i,j,k       Nw:Nw + 8,i,j,k
                            CCCFBAP(NwRe,i) = mean(CCCFBA(Nw:Nw + 4,i));
                            NwRe = NwRe + 1;
                        end
                    end
                %end
            %end
            %}
            %Real mean CCR (Armonics)
            %{
            [m1,n1,o1] = size(CCSFBA);
            %for k = 1:p1
                for j = 1:o1
                    %{
                    for i = 1:m1
                        CCSFBAW(i,j,k) = mean(CCSFBA(i,:,j,k));
                        CCCFBAW(i,j,k) = mean(CCCFBA(i,:,j,k));
                    end
                    %}
                    for i = 1:n1
                        NwRe = 1;
                        for Nw = 1:8:m1 - 4    %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
                            CCSFBAP(NwRe,i,j) = mean(CCSFBA(Nw:Nw + 4,i,j));   %Nw:Nw + 2,i,j,k      Nw:Nw + 4,i,j,k       Nw:Nw + 4,i,j,k       Nw:Nw + 8,i,j,k
                            CCCFBAP(NwRe,i,j) = mean(CCCFBA(Nw:Nw + 4,i,j));
                            NwRe = NwRe + 1;
                        end
                    end
                end
            %end
            %}
            %Real mean Filter Banck Windows (Armonics) CCRFBAP[NRealProbe,NArmonic,NFilterBank,NNomFrec]
            %
            [m1,n1,o1,p1] = size(CCSFBA);
            for k = 1:p1
                for j = 1:o1
                    %{
                    for i = 1:m1
                        CCSFBAW(i,j,k) = mean(CCSFBA(i,:,j,k));
                        CCCFBAW(i,j,k) = mean(CCCFBA(i,:,j,k));
                    end
                    %}
                    for i = 1:n1
                        NwRe = 1;
                        for Nw = 1:8:m1 - 4    %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
                            CCSFBAP(NwRe,i,j,k) = mean(CCSFBA(Nw:Nw + 4,i,j,k));   %Nw:Nw + 2,i,j,k      Nw:Nw + 4,i,j,k       Nw:Nw + 4,i,j,k       Nw:Nw + 8,i,j,k
                            CCCFBAP(NwRe,i,j,k) = mean(CCCFBA(Nw:Nw + 4,i,j,k));
                            NwRe = NwRe + 1;
                        end
                    end
                end
            end
            %
            %Weights per filter
            for i = 1:6
                W(1,i) = (i)^(-1.25) + 0.25;
            end
            %Matriz size
            %[m2,n2] = size(CCSFBAP);
            %[m2,n2,o2] = size(CCSFBAP);
            [m2,n2,o2,p2] = size(CCSFBAP);
            %Sum FBCCR (Armonics)   RokXX = [NRealProbe,NArmonic,NNomFrec]
            %
            for k = 1:o2
                RoNAS(:,:,k,:) = W(1,k).*CCSFBAP(:,:,k,:).^2;
                RoNAC(:,:,k,:) = W(1,k).*CCCFBAP(:,:,k,:).^2;
            end
            for k = 1:p2
                for j = 1:n2
                    for i = 1:m2
                        RokAS(i,j,k) = sum(RoNAS(i,j,:,k));
                        RokAC(i,j,k) = sum(RoNAC(i,j,:,k));
                        RokAP(i,j,k) = mean([RokAS(i,j,k) RokAC(i,j,k)]);
                    end
                end
            end
            %
            %Clasification CCR
            %{
            for i = 1:m2
                for j = 1:n2
                    CCSFBAP(i,j) = mean([CCSFBAP(i,j) CCCFBAP(i,j)]);
                end
                %for k = 1:o2
                    %[StimFBAS(i,k),ValueFBAS(i,k)] = MaxCan(RokAS(i,:,k),fxA(k,:));
                    %[StimFBAC(i,k),ValueFBAC(i,k)] = MaxCan(RokAC(i,:,k),fxA(k,:));
                    
                    %[StimFBAP(i,k),ValueFBAP(i,k)] = MaxCan(RokAP(i,:,k),fxA(k,:));
                %end
                [FinalStimS,FinalValueS] = MaxCan(CCSFBAP(i,:),fx);
                [FinalStimC,FinalValueC] = MaxCan(CCCFBAP(i,:),fx);
                
                [FinalStimP,FinalValueP] = MaxCan(CCSFBAP(i,:),fx);
                
                laCannonS.(Sub)(NS,ses,i) = FinalStimS;
                laCannonC.(Sub)(NS,ses,i) = FinalStimC;
                laCannonP.(Sub)(NS,ses,i) = FinalStimP;
            end
            %}
            %Clasification CCR (Armonics)
            %{
            for i = 1:m2
                for j = 1:n2
                    CCSFBAP(i,j) = mean([CCSFBAP(i,j) CCCFBAP(i,j)]);
                end
                for k = 1:o2
                    [StimFBAS(i,k),ValueFBAS(i,k)] = MaxCanA(CCSFBAP(i,:,k),fxA(k,:));
                    [StimFBAC(i,k),ValueFBAC(i,k)] = MaxCanA(CCCFBAP(i,:,k),fxA(k,:));
                    
                    [StimFBAP(i,k),ValueFBAP(i,k)] = MaxCanA(CCSFBAP(i,:,k),fxA(k,:));
                end
                [FinalStimS,FinalValueS] = MaxCan(ValueFBAS(i,:),StimFBAS(i,:));
                [FinalStimC,FinalValueC] = MaxCan(ValueFBAC(i,:),StimFBAC(i,:));
                
                [FinalStimP,FinalValueP] = MaxCan(ValueFBAP(i,:),StimFBAP(i,:));
                
                laCannonS.(Sub)(NS,ses,i) = FinalStimS;
                laCannonC.(Sub)(NS,ses,i) = FinalStimC;
                laCannonP.(Sub)(NS,ses,i) = FinalStimP;
            end
            %}
            %Clasification FBCCR (Armonics)
            %
            for i = 1:m2
                for k = 1:o2
                    [StimFBAS(i,k),ValueFBAS(i,k)] = MaxCanA(RokAS(i,:,k),fxA(k,:));
                    [StimFBAC(i,k),ValueFBAC(i,k)] = MaxCanA(RokAC(i,:,k),fxA(k,:));
                    
                    [StimFBAP(i,k),ValueFBAP(i,k)] = MaxCanA(RokAP(i,:,k),fxA(k,:));
                end
                [FinalStimS,FinalValueS] = MaxCan(ValueFBAS(i,:),StimFBAS(i,:));
                [FinalStimC,FinalValueC] = MaxCan(ValueFBAC(i,:),StimFBAC(i,:));
                
                [FinalStimP,FinalValueP] = MaxCan(ValueFBAP(i,:),StimFBAP(i,:));
                
                laCannonS.(Sub)(NS,ses,i) = FinalStimS;
                laCannonC.(Sub)(NS,ses,i) = FinalStimC;
                laCannonP.(Sub)(NS,ses,i) = FinalStimP;
            end
            %
            
            %-------------------------------------------------------------%
            %DSPv3
            %Xx = Data.TargetFrequency(NS);
            
            %EEG = Data.EEG(:,NS);
            %fn = Data.AmpSamlingFrequency;                  %Frecuencia de muestreo
            %Ts = 1/fn;                                      %Perdiodo de muestreo
            N = length(EEG);                                % # de muestras de la señal
            
            %f = (0:N-1)*fn/N;                              % Espectro de la ventana
            
            for windowss = 1:m1
                %SWW2=Wds.';
                Tfour = abs(fft(Wds(windowss,:)'));
                if rem(N,2) ~= 0
                    selec = 1 : (N+1)/2;    %Si es impar
                else
                    selec = 1 : N/2 + 1;    %Si es Par
                end
                
                %Densidad espectar la potencia
                
                DSP =(Tfour.^2)/(N*fn);
                DSPg = DSP(selec);
                DSP_g = [DSPg(1); 2*DSPg(2:end-1); DSPg(end)];
                f_g = f(selec);
                
                [Res, Fff] = max(DSP_g);
                R = Fff * 0.062;                                   %Frecuencia
                
                %if (Xx-2.5) < R && R < (Xx+2.5)
                    %R = Xx;
                %end
                
                Promedio(1,windowss) = R;
                
            end
            NwRe = 1;
            for Nw = 1:8:m1 - 4       %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
                [elmejor, dig]= Valores_repetidos(Promedio(1,Nw:Nw + 4));   %1,Nw:Nw + 2    1,Nw:Nw + 4     1,Nw:Nw + 4     1,Nw:Nw + 8
                laDig.(Sub)(NS,ses,NwRe) = dig;
                NwRe = NwRe + 1;
            end
        end
    end
end
toc

%% Stadistics
clearvars -except laCannonS laCannoC laCannonP laDig data
frecTest = [7.5,8.2,7,8.2,6,7.5,6,6,8.2,8.2,...
    9.3,6,6,8.2,6,6.5,7.5,7,6,6];
Yy = vertcat(frecTest', frecTest', frecTest', frecTest', frecTest', frecTest', frecTest', frecTest');
for sub = 1:5
    Sub = strcat("Sub",num2str(sub));
    allCannon(:,sub) = vertcat(laCannonP.(Sub)(:,1,1), laCannonP.(Sub)(:,2,1), ...
        laCannonP.(Sub)(:,1,2), laCannonP.(Sub)(:,2,2), ...
        laCannonP.(Sub)(:,1,3), laCannonP.(Sub)(:,2,3), ...
        laCannonP.(Sub)(:,1,4), laCannonP.(Sub)(:,2,4), ...
        laCannonP.(Sub)(:,1,5), laCannonP.(Sub)(:,2,5), ...
        laCannonP.(Sub)(:,1,6), laCannonP.(Sub)(:,2,6), ...
        laCannonP.(Sub)(:,1,7), laCannonP.(Sub)(:,2,7), ...
        laCannonP.(Sub)(:,1,8), laCannonP.(Sub)(:,2,8));
    allPSD(:,sub) = vertcat(laDig.(Sub)(:,1,1), laDig.(Sub)(:,2,1), ...
        laDig.(Sub)(:,1,2), laDig.(Sub)(:,2,2), ...
        laDig.(Sub)(:,1,3), laDig.(Sub)(:,2,3), ...
        laDig.(Sub)(:,1,4), laDig.(Sub)(:,2,4), ...
        laDig.(Sub)(:,1,5), laDig.(Sub)(:,2,5), ...
        laDig.(Sub)(:,1,6), laDig.(Sub)(:,2,6), ...
        laDig.(Sub)(:,1,7), laDig.(Sub)(:,2,7), ...
        laDig.(Sub)(:,1,8), laDig.(Sub)(:,2,8));
    for i = 1:160
        if (allCannon(i,sub) == Yy(i,1)) || (allCannon(i,sub) == Yy(i,1)*2) || (allCannon(i,sub) == Yy(i,1)*3) || (allCannon(i,sub) == Yy(i,1)*4) || (allCannon(i,sub) == Yy(i,1)*5)
            contCCR(1,i) = 1;
        else
            contCCR(1,i) = 0;
        end
        if ((allPSD(i,sub) <= Yy(i,1) + 0.25) && (allPSD(i,sub) >= Yy(i,1) - 0.25)) || (allPSD(i,sub) <= ((Yy(i,1)*2) + 0.25) && (allPSD(i,sub) >= (Yy(i,1)*2) - 0.25)) || ((allPSD(i,sub) <= (Yy(i,1)*3) + 0.25) && (allPSD(i,sub) >= (Yy(i,1)*3) - 0.25)) || ((allPSD(i,sub) <= (Yy(i,1)*4) + 0.25) && (allPSD(i,sub) >= (Yy(i,1)*4) - 0.25)) || ((allPSD(i,sub) <= (Yy(i,1)*5) + 0.25) && (allPSD(i,sub) >= (Yy(i,1)*5) - 0.25))
            contPSD(1,i) = 1;
        else
            contPSD(1,i) = 0;
        end
    end
    ExactCCR(1,sub) = sum(contCCR)/160;
    ExactPSD(1,sub) = sum(contPSD)/160;
end
