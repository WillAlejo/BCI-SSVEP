%% Characteristics
clear;
clc;

tic
%for sub = 1:5
    
    sub = 5;
    %for ses = 1:2
        
        ses = 1;
        %Load Database (specific values)
        load("Sub" + sub + "_" + ses + "_multitarget.mat")
        
        fn = Data.AmpSamlingFrequency;
        Sub = strcat("Sub",num2str(sub));
        %for NS = 1:10
            
            NS = 1;
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
                laCannonP.(Sub)(NS,ses,i) = FinalValueP;
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
        %end
    %end
%end
toc
%% Data Training
clearvars -except laCannonS laCannoC laCannonP laDig data
%
allCannon = vertcat(laCannonP.Sub1(:,1,1),laCannonP.Sub1(:,2,1),...
    laCannonP.Sub1(:,1,2),laCannonP.Sub1(:,2,2),...
    laCannonP.Sub1(:,1,3),laCannonP.Sub1(:,2,3),...
    laCannonP.Sub1(:,1,4),laCannonP.Sub1(:,2,4),...
    laCannonP.Sub1(:,1,5),laCannonP.Sub1(:,2,5),...
    laCannonP.Sub1(:,1,6),laCannonP.Sub1(:,2,6),...
    laCannonP.Sub1(:,1,7),laCannonP.Sub1(:,2,7),...
    laCannonP.Sub1(:,1,8),laCannonP.Sub1(:,2,8),...
    laCannonP.Sub2(:,1,1),laCannonP.Sub2(:,2,1),...
    laCannonP.Sub2(:,1,2),laCannonP.Sub2(:,2,2),...
    laCannonP.Sub2(:,1,3),laCannonP.Sub2(:,2,3),...
    laCannonP.Sub2(:,1,4),laCannonP.Sub2(:,2,4),...
    laCannonP.Sub2(:,1,5),laCannonP.Sub2(:,2,5),...
    laCannonP.Sub2(:,1,6),laCannonP.Sub2(:,2,6),...
    laCannonP.Sub2(:,1,7),laCannonP.Sub2(:,2,7),...
    laCannonP.Sub2(:,1,8),laCannonP.Sub2(:,2,8),...
    laCannonP.Sub3(:,1,1),laCannonP.Sub3(:,2,1),...
    laCannonP.Sub3(:,1,2),laCannonP.Sub3(:,2,2),...
    laCannonP.Sub3(:,1,3),laCannonP.Sub3(:,2,3),...
    laCannonP.Sub3(:,1,4),laCannonP.Sub3(:,2,4),...
    laCannonP.Sub3(:,1,5),laCannonP.Sub3(:,2,5),...
    laCannonP.Sub3(:,1,6),laCannonP.Sub3(:,2,6),...
    laCannonP.Sub3(:,1,7),laCannonP.Sub3(:,2,7),...
    laCannonP.Sub3(:,1,8),laCannonP.Sub3(:,2,8),...
    laCannonP.Sub4(:,1,1),laCannonP.Sub4(:,2,1),...
    laCannonP.Sub4(:,1,2),laCannonP.Sub4(:,2,2),...
    laCannonP.Sub4(:,1,3),laCannonP.Sub4(:,2,3),...
    laCannonP.Sub4(:,1,4),laCannonP.Sub4(:,2,4),...
    laCannonP.Sub4(:,1,5),laCannonP.Sub4(:,2,5),...
    laCannonP.Sub4(:,1,6),laCannonP.Sub4(:,2,6),...
    laCannonP.Sub4(:,1,7),laCannonP.Sub4(:,2,7),...
    laCannonP.Sub4(:,1,8),laCannonP.Sub4(:,2,8),...
    laCannonP.Sub5(:,1,1),laCannonP.Sub5(:,2,1),...
    laCannonP.Sub5(:,1,2),laCannonP.Sub5(:,2,2),...
    laCannonP.Sub5(:,1,3),laCannonP.Sub5(:,2,3),...
    laCannonP.Sub5(:,1,4),laCannonP.Sub5(:,2,4),...
    laCannonP.Sub5(:,1,5),laCannonP.Sub5(:,2,5),...
    laCannonP.Sub5(:,1,6),laCannonP.Sub5(:,2,6),...
    laCannonP.Sub5(:,1,7),laCannonP.Sub5(:,2,7),...
    laCannonP.Sub5(:,1,8),laCannonP.Sub5(:,2,8));
%
%{
allCannon = vertcat(laCannonP.Sub1(:,1,1),laCannonP.Sub1(:,2,1),...
    laCannonP.Sub1(:,1,2),laCannonP.Sub1(:,2,2),...
    laCannonP.Sub1(:,1,3),laCannonP.Sub1(:,2,3),...
    laCannonP.Sub1(:,1,4),laCannonP.Sub1(:,2,4),...
    laCannonP.Sub2(:,1,1),laCannonP.Sub2(:,2,1),...
    laCannonP.Sub2(:,1,2),laCannonP.Sub2(:,2,2),...
    laCannonP.Sub2(:,1,3),laCannonP.Sub2(:,2,3),...
    laCannonP.Sub2(:,1,4),laCannonP.Sub2(:,2,4),...
    laCannonP.Sub3(:,1,1),laCannonP.Sub3(:,2,1),...
    laCannonP.Sub3(:,1,2),laCannonP.Sub3(:,2,2),...
    laCannonP.Sub3(:,1,3),laCannonP.Sub3(:,2,3),...
    laCannonP.Sub3(:,1,4),laCannonP.Sub3(:,2,4),...
    laCannonP.Sub4(:,1,1),laCannonP.Sub4(:,2,1),...
    laCannonP.Sub4(:,1,2),laCannonP.Sub4(:,2,2),...
    laCannonP.Sub4(:,1,3),laCannonP.Sub4(:,2,3),...
    laCannonP.Sub4(:,1,4),laCannonP.Sub4(:,2,4),...
    laCannonP.Sub5(:,1,1),laCannonP.Sub5(:,2,1),...
    laCannonP.Sub5(:,1,2),laCannonP.Sub5(:,2,2),...
    laCannonP.Sub5(:,1,3),laCannonP.Sub5(:,2,3),...
    laCannonP.Sub5(:,1,4),laCannonP.Sub5(:,2,4));
%}
%
allDSP = vertcat(laDig.Sub1(:,1,1),laDig.Sub1(:,2,1),...
    laDig.Sub1(:,1,2),laDig.Sub1(:,2,2),...
    laDig.Sub1(:,1,3),laDig.Sub1(:,2,3),...
    laDig.Sub1(:,1,4),laDig.Sub1(:,2,4),...
    laDig.Sub1(:,1,5),laDig.Sub1(:,2,5),...
    laDig.Sub1(:,1,6),laDig.Sub1(:,2,6),...
    laDig.Sub1(:,1,7),laDig.Sub1(:,2,7),...
    laDig.Sub1(:,1,8),laDig.Sub1(:,2,8),...
    laDig.Sub2(:,1,1),laDig.Sub2(:,2,1),...
    laDig.Sub2(:,1,2),laDig.Sub2(:,2,2),...
    laDig.Sub2(:,1,3),laDig.Sub2(:,2,3),...
    laDig.Sub2(:,1,4),laDig.Sub2(:,2,4),...
    laDig.Sub2(:,1,5),laDig.Sub2(:,2,5),...
    laDig.Sub2(:,1,6),laDig.Sub2(:,2,6),...
    laDig.Sub2(:,1,7),laDig.Sub2(:,2,7),...
    laDig.Sub2(:,1,8),laDig.Sub2(:,2,8),...
    laDig.Sub3(:,1,1),laDig.Sub3(:,2,1),...
    laDig.Sub3(:,1,2),laDig.Sub3(:,2,2),...
    laDig.Sub3(:,1,3),laDig.Sub3(:,2,3),...
    laDig.Sub3(:,1,4),laDig.Sub3(:,2,4),...
    laDig.Sub3(:,1,5),laDig.Sub3(:,2,5),...
    laDig.Sub3(:,1,6),laDig.Sub3(:,2,6),...
    laDig.Sub3(:,1,7),laDig.Sub3(:,2,7),...
    laDig.Sub3(:,1,8),laDig.Sub3(:,2,8),...
    laDig.Sub4(:,1,1),laDig.Sub4(:,2,1),...
    laDig.Sub4(:,1,2),laDig.Sub4(:,2,2),...
    laDig.Sub4(:,1,3),laDig.Sub4(:,2,3),...
    laDig.Sub4(:,1,4),laDig.Sub4(:,2,4),...
    laDig.Sub4(:,1,5),laDig.Sub4(:,2,5),...
    laDig.Sub4(:,1,6),laDig.Sub4(:,2,6),...
    laDig.Sub4(:,1,7),laDig.Sub4(:,2,7),...
    laDig.Sub4(:,1,8),laDig.Sub4(:,2,8),...
    laDig.Sub5(:,1,1),laDig.Sub5(:,2,1),...
    laDig.Sub5(:,1,2),laDig.Sub5(:,2,2),...
    laDig.Sub5(:,1,3),laDig.Sub5(:,2,3),...
    laDig.Sub5(:,1,4),laDig.Sub5(:,2,4),...
    laDig.Sub5(:,1,5),laDig.Sub5(:,2,5),...
    laDig.Sub5(:,1,6),laDig.Sub5(:,2,6),...
    laDig.Sub5(:,1,7),laDig.Sub5(:,2,7),...
    laDig.Sub5(:,1,8),laDig.Sub5(:,2,8));
%
%{
allDSP = vertcat(laDig.Sub1(:,1,1),laDig.Sub1(:,2,1),...
    laDig.Sub1(:,1,2),laDig.Sub1(:,2,2),...
    laDig.Sub1(:,1,3),laDig.Sub1(:,2,3),...
    laDig.Sub1(:,1,4),laDig.Sub1(:,2,4),...
    laDig.Sub2(:,1,1),laDig.Sub2(:,2,1),...
    laDig.Sub2(:,1,2),laDig.Sub2(:,2,2),...
    laDig.Sub2(:,1,3),laDig.Sub2(:,2,3),...
    laDig.Sub2(:,1,4),laDig.Sub2(:,2,4),...
    laDig.Sub3(:,1,1),laDig.Sub3(:,2,1),...
    laDig.Sub3(:,1,2),laDig.Sub3(:,2,2),...
    laDig.Sub3(:,1,3),laDig.Sub3(:,2,3),...
    laDig.Sub3(:,1,4),laDig.Sub3(:,2,4),...
    laDig.Sub4(:,1,1),laDig.Sub4(:,2,1),...
    laDig.Sub4(:,1,2),laDig.Sub4(:,2,2),...
    laDig.Sub4(:,1,3),laDig.Sub4(:,2,3),...
    laDig.Sub4(:,1,4),laDig.Sub4(:,2,4),...
    laDig.Sub5(:,1,1),laDig.Sub5(:,2,1),...
    laDig.Sub5(:,1,2),laDig.Sub5(:,2,2),...
    laDig.Sub5(:,1,3),laDig.Sub5(:,2,3),...
    laDig.Sub5(:,1,4),laDig.Sub5(:,2,4));
%}
[allCannon,cannonMin,cannonMax] = normalizacion(allCannon);
[allDSP,dspMin,dspMax] = normalizacion(allDSP);
allCharactDatasH = horzcat(allCannon,allDSP);

%General Classification
Resps = [7.5000 8.2000 7 8.2000 6 7.5000 6 6 8.2000 8.2000 9.3000 6 6 8.2000 6 6.5000 7.5000 7 6 6];
allTrainingDatas = allCharactDatasH;
CordX = allTrainingDatas;
%{
CordY = vertcat(Resps',Resps',Resps',Resps',Resps',Resps',Resps',Resps',...
    Resps',Resps',Resps',Resps',Resps',Resps',Resps',Resps',...
    Resps',Resps',Resps',Resps',Resps',Resps',Resps',Resps',...
    Resps',Resps',Resps',Resps',Resps',Resps',Resps',Resps',...
    Resps',Resps',Resps',Resps',Resps',Resps',Resps',Resps');
%}
%
CordY = vertcat(Resps',Resps',Resps',Resps',...
    Resps',Resps',Resps',Resps',...
    Resps',Resps',Resps',Resps',...
    Resps',Resps',Resps',Resps',...
    Resps',Resps',Resps',Resps');
%
plot(CordX(:,1),CordX(:,2),'x')

K = 1;
d = "euclidean";

%Md1 = fitcknn(CordX,CordY,'NumNeighbors',K,'Distance',d);
Md1 = fitcauto(CordX,CordY);
%% Save Classifier
STR.fullData = allTrainingDatas;
STR.bestModel = Md1;
STR.cannonMin = cannonMin;
STR.cannonMax = cannonMax;
STR.dspMin = dspMin;
STR.dspMax = dspMax;
save('Best_9W4S_Classifier','-struct','STR');