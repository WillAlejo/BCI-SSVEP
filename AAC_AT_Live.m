%% Data Adquisition
clear;
clc;
warning('off')

Ndatas=2.5;
maxduration = 1000;
%deadline = GetSecs + maxduration;

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
%while (~KbCheck) && (GetSecs < deadline)
for li = 1:10
    tic
    board_shim.start_stream(45000, '');
    %WaitSecs(Ndatas);
    pause(Ndatas)
    board_shim.stop_stream();
    %nfft = DataFilter.get_nearest_power_of_two(sampling_rate);
    data = board_shim.get_board_data(board_shim.get_board_data_count());  %Dos
    %data = board_shim.get_current_board_data(10);                        %Uno
    %disp(data);
    %function resistance_channels = get_resistance_channels(board_id)    %Tres

    %% Characteristics Adquisitation; Live Data
    %clearvars -except a board_descr board_id board_shim params Ndatas data
    %One Test
    %
    %try
        fn = double(board_descr.sampling_rate);
        EEGFG = data(2,2:fn*(Ndatas-0.5) + 1);
        EEG = EEGFG';

        t = 0 : 1/fn : ((length(EEG)) - 1)/fn;
        f = 0 : fn/length(EEG) : fn - (fn/length(EEG));

        EEGNT = NoTch(EEG,3,59,61,fn);
        EEGF = BandPass(EEGNT,6,5,50,fn);

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
        for i = 1:3
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
        for k = 1:p2
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
    %catch
%         fprintf('ERROR \n');
%         board_shim.release_session();
    %end
    %end
    %
    %% Classification
    %clearvars -except sampling_rate Ndatas data laCannonS laCannonC laCannonP laDig
    load('Best_5W2S_Classifier.mat')
    laCannonP=laCannonP'-cannonMin;
    laCannonP=laCannonP/cannonMax;
    laDig=laDig'-dspMin;
    laDig=laDig/dspMax;
    allTestDatas = horzcat(laCannonP,laDig);
    
    Predecir1 = predict(bestModel,allTestDatas);
    NET.addAssembly('System.Speech');
    obj = System.Speech.Synthesis.SpeechSynthesizer;
    obj.Volume = 100;
    
    xlk = num2str(Predecir1);
    
    Speak(obj, xlk);

    disp(xlk);

    %WaitSecs(0.1);

    %clearvars -except a board_descr board_id board_shim params Ndatas data deadline
    toc
end
board_shim.release_session();
fprintf('Turn Off \n');
