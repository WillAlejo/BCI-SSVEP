%% Data Adquisition
clear;
clc;

% Ndatas = 2;
% BoardShim.set_log_file('brainflow.log');
% BoardShim.enable_dev_board_logger();
% 
% params = BrainFlowInputParams();
% board_shim = BoardShim(int32(BoardIDs.CYTON_BOARD), params);
% board_id = int32(BoardIDs.CYTON_BOARD);
% board_descr = BoardShim.get_board_descr(board_id);
% sampling_rate = int32(board_descr.sampling_rate);
% board_shim.prepare_session();
% board_shim.start_stream(45000, '');
% pause(Ndatas);
% board_shim.stop_stream();
% 
% nfft = DataFilter.get_nearest_power_of_two(sampling_rate);
% data = board_shim.get_board_data(board_shim.get_board_data_count());
% board_shim.release_session();

Ndatas=20.5;
maxduration = 1000;
deadline = GetSecs + maxduration;

fprintf('Turn on \n');
fprintf('Staring. Wait a moment... \n');
%BoardShim.set_log_file('brainflow.log');
BoardShim.enable_dev_board_logger();

params = BrainFlowInputParams();
board_shim = BoardShim(int32(BoardIDs.CYTON_BOARD), params);
%board_id = int32(1);
board_id = int32(BoardIDs.CYTON_BOARD);
board_descr = BoardShim.get_board_descr(board_id);
board_shim.prepare_session();
a = board_shim.config_board('~6');
board_shim.start_stream(45000, '');
WaitSecs(Ndatas);
board_shim.stop_stream();
%nfft = DataFilter.get_nearest_power_of_two(sampling_rate);
data = board_shim.get_board_data(board_shim.get_board_data_count());  %Dos
%data = board_shim.get_current_board_data(10);                        %Uno
%disp(data);
%function resistance_channels = get_resistance_channels(board_id)    %Tres

%plot(dummyvar)
%% Characteristics Adquisitation; Live Data
clearvars -except sampling_rate Ndatas data board_descr
%One Test
%{
tic
fn = double(sampling_rate);
EEGFG = data(2,1:fn*Ndatas);
EEG = EEGFG';

t = 0 : 1/fn : ((length(EEG)) - 1)/fn;
f = 0 : fn/length(EEG) : fn - (fn/length(EEG));

EEGF = BandPass(EEG,6,5,50,fn);

%Windowing signal (Function) Wds = [NWindow,EEGWindow]
Wds = WindowingLive(EEGF,2,0.25);

%Frecuencies from CannonCorr
fx = [6 6.5 7 7.5 8.2 9.3];

%Frecuencies from CannonCorr (armonics)
fxA = zeros(1,5);
for i = 1:6
    fxA(i,:) = [fx(i) fx(i)*2 fx(i)*3 fx(i)*4 fx(i)*5];
end

%Filter Bank (Function) EEGFB = [NWindow,EEGWindow,NFilterBank]
LF = [4 10 16 22 28 34];
HF = [48.5 48.5 48.5 48.5 48.5 48.5];

EEGFB = FilterBank(Wds,6,LF,HF,fn);

%FBCCR (Armonics) Cannon[Sen,cos] =
%{Sen[NWindow,NArmonic,NFilterBank,NNomFrec],Cos[NWindow,NArmonic,NFilterBank,NNomFrec]}
for j = 1:6
    for i = 1:6
        [CCSFBA(:,:,i,j),CCCFBA(:,:,i,j)] = Cannon(EEGFB(:,:,i),fxA(j,:),t);
    end
end

%Real mean Filter Banck Windows (Armonics) CCRFBAP[NArmonic,NFilterBank,NNomFrec]
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
            
            for Nw = 1:16:m1 - 8    %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
                CCSFBAP(i,j,k) = mean(CCSFBA(Nw:Nw + 8,i,j,k));   %Nw:Nw + 2,i,j,k      Nw:Nw + 4,i,j,k       Nw:Nw + 4,i,j,k       Nw:Nw + 8,i,j,k
                CCCFBAP(i,j,k) = mean(CCCFBA(Nw:Nw + 8,i,j,k));
                
            end
        end
    end
end

%Weights per filter
for i = 1:6
    W(1,i) = (i)^(-1.25) + 0.25;
end

%Sum FBCCR (Armonics)   RokXX = [NArmonic,NNomFrec]
[n2,o2,p2] = size(CCSFBAP);
for k = 1:o2
    RoNAS(:,k,:) = W(1,k).*CCSFBAP(:,k,:).^2;
    RoNAC(:,k,:) = W(1,k).*CCCFBAP(:,k,:).^2;
end
for k = 1:p2
    for j = 1:n2
        %for i = 1:m2
            RokAS(j,k) = sum(RoNAS(j,:,k));
            RokAC(j,k) = sum(RoNAC(j,:,k));
            RokAP(j,k) = mean([RokAS(j,k) RokAC(j,k)]);
        %end
    end
end
%RokAP(1,i,k) = mean([RokAS(1,i,k) RokAC(1,i,k)]);

%Clasification FBCCR (Armonics)
%for i = 1:m2
    for k = 1:o2
        [StimFBAS(1,k),ValueFBAS(1,k)] = MaxCan(RokAS(:,k),fxA(k,:));
        [StimFBAC(1,k),ValueFBAC(1,k)] = MaxCan(RokAC(:,k),fxA(k,:));

        [StimFBAP,ValueFBAP] = MaxCan(RokAP(:,k),fxA(k,:));
    end
    [FinalStimS,FinalValueS] = MaxCan(ValueFBAS(1,:),StimFBAS(1,:));
    [FinalStimC,FinalValueC] = MaxCan(ValueFBAC(1,:),StimFBAC(1,:));

    [FinalStimP,FinalValueP] = MaxCan(ValueFBAP(1,:),StimFBAP(1,:));

    laCannonS = FinalStimS;
    laCannonC = FinalStimC;
    laCannonP = FinalStimP;
%end

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
    %N = fn*Ndatas/2;
    %N = 200*4;
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
    R = (Fff-1) * 1/Ndatas;                                   %Frecuencia
    %R = Fff;
    %if (Xx-2.5) < R && R < (Xx+2.5)
    %   R = Xx;
    %end

    Promedio(1,windowss) = R;

end

%for Nw = 1:16:m1 - 8       %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
    [elmejor, dig]= Valores_repetidos(Promedio(1,:));   %1,Nw:Nw + 2    1,Nw:Nw + 4     1,Nw:Nw + 4     1,Nw:Nw + 8
    laDig = dig;
%end
%}
%N Tests
%
fn = double(board_descr.sampling_rate);
%fn = 250;
%Ndatas = 60;
%fn = double(board_descr.sampling_rate);
EEGFG = data(2,2:fn*(Ndatas-0.5) + 1);
EEG = EEGFG';
%EEG = EEGT(4000:12000-1,:);

t = 0 : 1/fn : ((length(EEG)) - 1)/fn;
f = 0 : fn/length(EEG) : fn - (fn/length(EEG));

EEGNT = NoTch(EEG,3,59,61,fn);
EEGF = BandPass(EEGNT,6,5,50,fn);

%Windowing signal (Function) Wds = [NWindow,EEGWindow]
Wds = WindowingLive(EEGF,1,0.25);

%Frecuencies from CannonCorr
fx = [6 6.5 7 7.5 8.2 9.3];

%Frecuencies from CannonCorr (armonics)
fxA = zeros(1,5);
for i = 1:6
    fxA(i,:) = [fx(i) fx(i)*2 fx(i)*3 fx(i)*4 fx(i)*5];
end

%Filter Bank (Function) EEGFB = [NWindow,EEGWindow,NFilterBank]
LF = [4 10 16 22 28 34];
HF = [48.5 48.5 48.5 48.5 48.5 48.5];

EEGFB = FilterBank(Wds,6,LF,HF,fn);

%FBCCR (Armonics) Cannon[Sen,cos] =
%{Sen[NWindow,NArmonic,NFilterBank,NNomFrec],Cos[NWindow,NArmonic,NFilterBank,NNomFrec]}
for j = 1:6
    for i = 1:6
        [CCSFBA(:,:,i,j),CCCFBA(:,:,i,j)] = Cannon(EEGFB(:,:,i),fxA(j,:),t);
    end
end

%Real mean Filter Banck Windows (Armonics) CCRFBAP[NRealProbe,NArmonic,NFilterBank,NNomFrec]
[m1,n1,o1,p1] = size(CCSFBA);
for k = 1:p1
    for j = 1:o1
        for i = 1:n1
            NwRe = 1;
            for Nw = 1:4:m1 - 2   %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
                CCSFBAP(NwRe,i,j,k) = mean(CCSFBA(Nw:Nw + 2,i,j,k));   %Nw:Nw + 2,i,j,k      Nw:Nw + 4,i,j,k       Nw:Nw + 4,i,j,k       Nw:Nw + 8,i,j,k
                CCCFBAP(NwRe,i,j,k) = mean(CCCFBA(Nw:Nw + 2,i,j,k ));
                NwRe = NwRe + 1;
            end
        end
    end
end

%Weights per filter
for i = 1:6
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
    for k = 1:o2
        [StimFBAS(i,k),ValueFBAS(i,k)] = MaxCan(RokAS(i,:,k),fxA(k,:));
        [StimFBAC(i,k),ValueFBAC(i,k)] = MaxCan(RokAC(i,:,k),fxA(k,:));

        [StimFBAP(i,k),ValueFBAP(i,k)] = MaxCan(RokAP(i,:,k),fxA(k,:));
    end
    [FinalStimS,FinalValueS] = MaxCan(ValueFBAS(i,:),StimFBAS(i,:));
    [FinalStimC,FinalValueC] = MaxCan(ValueFBAC(i,:),StimFBAC(i,:));

    [FinalStimP,FinalValueP] = MaxCan(ValueFBAP(i,:),StimFBAP(i,:));

    laCannonS(1,i) = FinalStimS;
    laCannonC(1,i) = FinalStimC;
    laCannonP(1,i) = FinalStimP;
end

%-------------------------------------------------------------%
%DSPv3

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
    R = (Fff-1) * 1/40;                                   %Frecuencia

    Promedio(1,windowss) = R;

end
NwRe = 1;
for Nw = 1:4:m1 - 2     %1:4:m1 - 2     1:8:m1 - 4      1:8:m1 - 4     1:16:m1 - 8
    [elmejor, dig]= Valores_repetidos(Promedio(1,Nw:Nw + 2));   %1,Nw:Nw + 2    1,Nw:Nw + 4     1,Nw:Nw + 4     1,Nw:Nw + 8
    laDig(1,NwRe) = dig;
    NwRe = NwRe + 1;
end
%
%% Classification
clearvars -except sampling_rate Ndatas data laCannonS laCannonC laCannonP laDig
load('Best_5W2S_Classifier.mat')
laCannonP=laCannonP'-cannonMin;
laCannonP=laCannonP/cannonMax;
laDig=laDig'-dspMin;
laDig=laDig/dspMax;
allTestDatas = horzcat(laCannonP,laDig);
Predecir1 = predict(bestModel,allTestDatas);
%% Stadistics
clearvars -except sampling_rate Ndatas data laCannonS laCannonC laCannonP laDig Predecir1
frecTest = 7.5;
Yy = [frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest];
%Yy = [frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest frecTest];
[estaPrecision,estaError,estaExactitud]=eficacia(Predecir1,Yy',frecTest);
%% Save Results
clearvars -except sampling_rate Ndatas data laCannonS laCannonC laCannonP laDig Predecir1 frecTest estaPrecision estaError estaExactitud
Stadistics.results = Predecir1;
Stadistics.frecTest = frecTest;
Stadistics.precision = estaPrecision;
Stadistics.error = estaError;
Stadistics.accuaracy = estaExactitud;
save("Stadistics_5W2S_" + frecTest + "_bestClass.mat",'-struct','Stadistics');