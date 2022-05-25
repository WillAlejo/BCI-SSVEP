%% Data Adquisition
clear;
clc;

Ndatas = 40.5;
maxduration = 1000;
deadline = GetSecs + maxduration;

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
board_shim.release_session();
%% DSP validation
%DSPv3
%N Tests
%
%fn = double(board_descr.sampling_rate);
fn = 250;
%Ndatas = 60;
%fn = double(board_descr.sampling_rate);
EEGFG = data(2,2:fn*(Ndatas-0.5) + 1);
EEG = EEGFG';
%EEG = data.EEG;
%EEG = EEGT(4000:12000-1,:);

t = 0 : 1/fn : ((length(EEG)) - 1)/fn;
f = 0 : fn/length(EEG) : fn - (fn/length(EEG));

EEGNT = NoTch(EEG,3,59,61,fn);
EEGF = BandPass(EEGNT,6,5,50,fn);

%Windowing signal (Function) Wds = [NWindow,EEGWindow]
Wds = WindowingLive(EEGF,1,0.25);



N = length(EEG);                                % # de muestras de la señal

%f = (0:N-1)*fn/N;                              % Espectro de la ventana

for windowss = 1:77
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
    DSP_g(:,windowss) = [DSPg(1); 2*DSPg(2:end-1); DSPg(end)];
    f_g = f(selec);

    [Res, Fff] = max(DSP_g(:,windowss));
    R = (Fff-1) * 1/40;                                   %Frecuencia

    Promedio(1,windowss) = R;
end
plot(f_g,DSP_g(:,50))
%% Data save
fn = double(board_descr.sampling_rate);
EEGFG = data(2,2:fn*(Ndatas-0.5) + 1);
EEG = EEGFG';
% EEG = 1;
% fn = 50;

subjet3.Data.EEG = EEG;
subjet3.Data.sampleRate = fn;

save('Subjet3_1_2','-struct','subjet3');