%% Characteristics Adquisitation
clear;
clc;
tic
for sub = 1:3

    %sub = 1;
    for ses = 1:4

        %ses = 1;
        %Load Database (specific values)
        load("Subjet" + sub + "_" + ses + ".mat")

        fn = Data.sampleRate;
        Sub = strcat("Sub",num2str(sub));

        %NS = 7;
        EEG = Data.EEG;

        %Time and frecuency vectors
        t = 0 : 1/fn : ((length(EEG)) - 1)/fn;
        f = 0 : fn/length(EEG) : fn - (fn/length(EEG));

        %FFT
        Sfft = ((abs(fft(EEG))).^2)./(length(EEG));

        %Band Pass Filter from database (Function)
        EEGNT = NoTch(EEG,3,59,61,fn);
        EEGF = BandPass(EEGNT,6,5,50,fn);

        %FFT filtered database
        SfftF = ((abs(fft(EEGF))).^2)./(length(EEG));

        %Windowing signal (Function) Wds = [NWindow,EEGWindow]
        Wds = WindowingLive(EEGF,2,0.25);

        %Frecuencies from CannonCorr
        fx = [6 7.2 8 9.6];

        %Frecuencies from CannonCorr (armonics)
        fxA = zeros(1,3);
        for i = 1:4
            fxA(i,:) = [fx(i) fx(i)*2 fx(i)*3];
        end

        %Filter Bank (Function) EEGFB = [NWindow,EEGWindow,NFilterBank]
        LF = [4 10 16 22 28 34];
        HF = [48.5 48.5 48.5 48.5 48.5 48.5];

        EEGFB = FilterBank(Wds,6,LF,HF,fn);

        %FBCCR (Armonics) Cannon[Sen,cos] =
        %{Sen[NWindow,NArmonic,NFilterBank,NNomFrec],Cos[NWindow,NArmonic,NFilterBank,NNomFrec]}
        for j = 1:4
            for i = 1:3
                [CCSFBA(:,:,i,j),CCCFBA(:,:,i,j)] = Cannon(EEGFB(:,:,i),fxA(j,:),t);
            end
        end

        %Real mean Filter Banck Windows (Armonics) CCRFBAP[NRealProbe,NArmonic,NFilterBank,NNomFrec]
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
                    for Nw = 1:16:m1 - 8    %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
                        CCSFBAP(NwRe,i,j,k) = mean(CCSFBA(Nw:Nw + 8,i,j,k));   %Nw:Nw + 2,i,j,k      Nw:Nw + 4,i,j,k       Nw:Nw + 4,i,j,k       Nw:Nw + 8,i,j,k
                        CCCFBAP(NwRe,i,j,k) = mean(CCCFBA(Nw:Nw + 8,i,j,k));
                        NwRe = NwRe + 1;
                    end
                end
            end
        end

        %Weights per filter
        for i = 1:3
            W(1,i) = (i)^(-1.25) + 0.25;
        end

        %Sum FBCCR (Armonics)   RokXX = [NRealProbe,NArmonic,NNomFrec]
        [m2,n2,o2,p2] = size(CCSFBAP);
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
        %RokAP(1,i,k) = mean([RokAS(1,i,k) RokAC(1,i,k)]);

        %Clasification FBCCR (Armonics)
        for i = 1:m2
            for k = 1:p2
                [StimFBAS(i,k),ValueFBAS(i,k)] = MaxCan(RokAS(i,:,k),fxA(k,:));
                [StimFBAC(i,k),ValueFBAC(i,k)] = MaxCan(RokAC(i,:,k),fxA(k,:));

                [StimFBAP(i,k),ValueFBAP(i,k)] = MaxCan(RokAP(i,:,k),fxA(k,:));
            end
            [FinalStimS,FinalValueS] = MaxCan(ValueFBAS(i,:),StimFBAS(i,:));
            [FinalStimC,FinalValueC] = MaxCan(ValueFBAC(i,:),StimFBAC(i,:));

            [FinalStimP,FinalValueP] = MaxCan(ValueFBAP(i,:),StimFBAP(i,:));

            laCannonS.(Sub)(ses,i) = FinalStimS;
            laCannonC.(Sub)(ses,i) = FinalStimC;
            laCannonP.(Sub)(ses,i) = FinalStimP;
        end
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
        for Nw = 1:16:m1 - 8       %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
            [elmejor, dig]= Valores_repetidos(Promedio(1,Nw:Nw + 8));   %1,Nw:Nw + 2    1,Nw:Nw + 4     1,Nw:Nw + 4     1,Nw:Nw + 8
            laDig.(Sub)(ses,NwRe) = dig;
            NwRe = NwRe + 1;
        end
    end
end
toc
%% Data Training
clearvars -except laCannonS laCannoC laCannonP laDig data

% allCannon = vertcat(laCannonP.Sub1(1,:)',laCannonP.Sub2(1,:)',laCannonP.Sub3(1,:)',...
%     laCannonP.Sub1(2,:)',laCannonP.Sub2(2,:)',laCannonP.Sub3(2,:)',...
%     laCannonP.Sub1(3,:)',laCannonP.Sub2(3,:)',laCannonP.Sub3(3,:)',...
%     laCannonP.Sub1(4,:)',laCannonP.Sub2(4,:)',laCannonP.Sub3(4,:)');
% 
% allDSP = vertcat(laDig.Sub1(1,:)',laDig.Sub2(1,:)',laDig.Sub3(1,:)',...
%     laDig.Sub1(2,:)',laDig.Sub2(2,:)',laDig.Sub3(2,:)',...
%     laDig.Sub1(3,:)',laDig.Sub2(3,:)',laDig.Sub3(3,:)',...
%     laDig.Sub1(4,:)',laDig.Sub2(4,:)',laDig.Sub3(4,:)');

allCannon = vertcat(laCannonP.Sub1(1,:)',laCannonP.Sub1(2,:)',laCannonP.Sub1(3,:)',laCannonP.Sub1(4,:)');
allDSP = vertcat(laDig.Sub1(1,:)',laDig.Sub1(2,:)',laDig.Sub1(3,:)',laDig.Sub1(4,:)');

[allCannon,cannonMin,cannonMax] = normalizacion(allCannon);
[allDSP,dspMin,dspMax] = normalizacion(allDSP);
allCharactDatasH = horzcat(allCannon,allDSP);

%General Classification
%{
Resps = [6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 ... 
    7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 ...
    8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 ...
    9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6];
%}
% %
% Resps = [6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 ... 
%     7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 ...
%     8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 ...
%     9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6];
% %
Resps = [6 6 6 6 6 6 6 6 6 6 ... 
     7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2 ...
     8 8 8 8 8 8 8 8 8 8 ...
     9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6 9.6];
allTrainingDatas = allCharactDatasH;
CordX = allTrainingDatas;

CordY = vertcat(Resps');
%
%plot(CordX(1:30,1),CordX(1:30,2),'x', CordX(31:60,1),CordX(31:60,2),'o', CordX(61:90,1),CordX(61:90,2),'O', CordX(91:120,1),CordX(91:120,2),'X')

%Md1 = fitcknn(CordX,CordY,'NumNeighbors',K,'Distance',d);
Md1 = fitcauto(CordX,CordY);
%% Save Classifier
STR.fullData = allTrainingDatas;
STR.fullTarjets = CordY;
STR.bestModel = Md1;
STR.cannonMin = cannonMin;
STR.cannonMax = cannonMax;
STR.dspMin = dspMin;
STR.dspMax = dspMax;
save('Best_9W4S_Classifier_Live_4','-struct','STR');