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
%% Data Training
clearvars -except laCannonS laCannoC laCannonP laDig data
frecTest = [7.5,8.2,7,8.2,6,7.5,6,6,8.2,8.2,...
    9.3,6,6,8.2,6,6.5,7.5,7,6,6];
noClass = ["NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC" "NC"];
Yy = vertcat(frecTest', frecTest', frecTest', frecTest', frecTest', frecTest', frecTest', frecTest');
NC = vertcat(noClass', noClass', noClass', noClass', noClass', noClass', noClass', noClass');
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
    allPSD(:,sub) = vertcat(laDig.(Sub)(:,2,1), laDig.(Sub)(:,1,1), ...
        laDig.(Sub)(:,1,2), laDig.(Sub)(:,2,2), ...
        laDig.(Sub)(:,1,3), laDig.(Sub)(:,2,3), ...
        laDig.(Sub)(:,1,4), laDig.(Sub)(:,2,4), ...
        laDig.(Sub)(:,1,5), laDig.(Sub)(:,2,5), ...
        laDig.(Sub)(:,1,6), laDig.(Sub)(:,2,6), ...
        laDig.(Sub)(:,1,7), laDig.(Sub)(:,2,7), ...
        laDig.(Sub)(:,1,8), laDig.(Sub)(:,2,8));
end
allCannon = vertcat(allCannon(:,1), allCannon(:,4), allCannon(:,3), allCannon(:,4), allCannon(:,5));
allPSD = vertcat(allPSD(:,1), allPSD(:,2), allPSD(:,3), allPSD(:,4), allPSD(:,5));

[allCannon,cannonMin,cannonMax] = normalizacion(allCannon);
[allPSD,dspMin,dspMax] = normalizacion(allPSD);

%allTrainingDatas = horzcat(allCannon(1:640,:), allPSD(1:640,:));
%allTestDatas = horzcat(allCannon(641:800,:), allPSD(641:800,:));
allTrainingDatas = horzcat(allCannon, allPSD);

CordX = allTrainingDatas;
%CordY = [Yy; Yy; Yy; Yy];
CordY = [Yy; Yy; Yy; Yy; NC];
%CordY = [Yy; Yy; Yy; Yy;Yy];
%plot(CordX(:,1),CordX(:,2),'x',allTestDatas(:,1),allTestDatas(:,2),'O')
%plot(CordX(:,1),CordX(:,2),'x')
gscatter(CordX(:,1),CordX(:,2),CordY,'kmcrgbk','o*sd^hx',[],'on','FBCCA','PSD')

%Md1 = fitcauto(CordX,CordY);
%% Classifier validation
Predecir1 = predict(Md1,allTestDatas);

[prec1,error1,exactitud1]=eficacia(Predecir1,Yy,6);
[prec2,error2,exactitud2]=eficacia(Predecir1,Yy,6.5);
[prec3,error3,exactitud3]=eficacia(Predecir1,Yy,7);
[prec4,error4,exactitud4]=eficacia(Predecir1,Yy,7.5);
[prec5,error5,exactitud5]=eficacia(Predecir1,Yy,8.2);
[prec6,error6,exactitud6]=eficacia(Predecir1,Yy,9.3);

finalPrec = (prec1 + prec2 + prec3 + prec4 + prec5 + prec6)/6;
finalError = (error1 + error2 + error3 + error4 + error5 + error6)/6;
finalExactitud = (exactitud1 + exactitud2 + exactitud3 + exactitud4 + exactitud5 + exactitud6)/6;
%% Save Classifier
STR.fullData = allTrainingDatas;
STR.fullTarjets = CordY;
STR.bestModel = Md1;
STR.cannonMin = cannonMin;
STR.cannonMax = cannonMax;
STR.dspMin = dspMin;
STR.dspMax = dspMax;
save('Best_5W2S_Classifier_2.mat','-struct','STR');