function Estimulos() 
clear; clc; close all;
%clearvars -except sampling_rate Ndatas data
%%Initiate frequency
% Frames Period Freq. Simulated signal. 0 light. 1 dark
% [#]
% [ms] [Hz] [-]
% 3.0 50 .00 20.00 011
% 4.0 66.67 15.00 0011
% 5.0 83.33 12.00 00111
% 6.0 100.00 10.00 000111 
% 7.0 116.67 8.57 0001111
% 8.0 133.33 7.50 00001111
% 9.0 150.00 6.66 0000 11111

%According to the paper 1 is blackbox, 0 is white

% Seis =           [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];  %6          
% SeisCinco =      [0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1];      %6.5 
% SieteCinco =     [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1];            %7.5
% Doce =           [0 0 0 0 0 0 1 1 1 1 1 1];                          %12 
% 
% freq{1} = SieteCinco;   %4/8                 Arriba
% freq{2} = Seis;    %3/6                      Derecha
% freq{3} = Doce; %3/5                         Abajo
% freq{4} = SeisCinco; %4/7                    Izquierda




%9.6, 8, 6.8, y 6 :v

Nueve_Seis =     [0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];                    %9.6
Ocho =           [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1];              %8
Siete_Dos =      [0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1];          %7.2
Seis =           [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];  %6

freq{1} = Nueve_Seis;   %4/8                Arriba
freq{2} = Ocho;    %3/6                     Derecha
freq{3} = Siete_Dos; %3/5                   Abajo
freq{4} = Seis; %4/7                        Izquierda





% initiate freq table



%zzzzz= size(freq,2);

%%Generate display matrixes for movies
% Find LCM of freq matrix to create equal matrixes for all freqs 
lcmFreq = lcms([length(freq{1}),length(freq{2}),length(freq{3}),length(freq{4})]);
%Generate full movie matrix of frequency 
for i=1:4
freqCombine(i,:) = repmat(freq{i},1,lcmFreq/length(freq{i})); 
end
%Revert value because in Matlab 255 is white and 0 is black
freqCombine = 1 - freqCombine;

maxduration = 1000;


if nargin < 2
    frequency = 5;
end



% indexflip =1;
%for sss = 1:5
%    indexflip =indexflip+1;
%    textureValue = freqCombine(:, indexflip) .* [1; 2; 4; 8];
%    textureValue = textureValue(4)+textureValue(3)+textureValue(2)+ textureValue(1) +1;
%end
% cscsc =1;
 
 
 
 
try
    Screen('Preference','SkipSyncTests',1)
    myScreen = max(Screen('Screens'));
    %Priority(0);
    %myScreen = 2;
    [win,winRect] =   Screen('OpenWindow',myScreen,[]);
        
    [width, height] = RectSize(winRect);
    
    % Background color dark green, just to make sure
    Screen('FillRect',win,[0 127 0]);
    %%Make movie 
    targetWidth = 250;
    targetHeight = 250;
    
    % make textures clipped to screen size
    % Draw texture to screen: Draw 16 states or texture depens on the value of
    screenMatrix = flickerTexture(width, height, targetWidth, targetHeight);

    MatrizOscura = ones(height,width,'uint8');
    MatrizOscuraT = Screen('MakeTexture', win, MatrizOscura*255);

    for  i = 1:16 
    texture(i) = Screen('MakeTexture', win, uint8(screenMatrix{i})*255);
    end
    
   
    % Define refresh rate.
    ifi = Screen('GetFlipInterval', win);
    
   %Run in this duration
   
    deadline = GetSecs + maxduration;
        
    % Preview texture briefly before flickering
    % n.b. here we draw to back buffer
    Screen('DrawTexture',win,texture(16)); 
    VBLTimestamp = Screen('Flip', win, ifi);
    WaitSecs(2);
    
    % loop swapping buffers, checking keyboard, and checking time
    % param 2 denotes "dont clear buffer on flip", i.e., we alternate
    % our buffers cum textures
    indexflip = 1;
    % textureValue =0;    
    halfifi = 0.5*ifi;
    vbl =0; 
    C1 = 0;
    
%% Start looping movie   
    Priority(1);
    while (~KbCheck)&&(GetSecs < deadline)



    % Drawing
    %Compute texture value based on display value from freq long matrixes
    textureValue = freqCombine(:, indexflip) .* [1; 2; 4; 8];   
    textureValue = textureValue(4)+textureValue(3)+textureValue(2)+ textureValue(1) +1;
    %Draw it on the back buffer
    Screen('DrawTexture',win,texture(textureValue)); 
    %Display current index
    %Screen('DrawText', win, num2str(indexflip),400,400, 255);
    %Tell PTB no more drawing commands will be issued until the next flip
    
    % Fliping     
    %Screen('Flip', win, vbl + halfifi);
    %Flip ASAP
    [nx, ny, bbox] = DrawFormattedText(win,'¿Cómo te encuentras?',5,540,[255 0 0], 30);             %Izquierda
    [nx, ny, bbox] = DrawFormattedText(win, 'Tengo hambre', 880, 125, [255 0 0], 30);               %Arriba
    [nx, ny, bbox] = DrawFormattedText(win, '¡Hola!', 1765, 540, [255 0 0], 30);                    %Derecha
    [nx, ny, bbox] = DrawFormattedText(win, '¡Quiero ir al baño!', 860, 955, [255 0 0], 30);        %Abajo
    
   
    Screen('DrawingFinished', win);
    Screen('Flip', win);
    
    indexflip = indexflip+1;
    


    %Reset index at the end of freq matrix
        if indexflip > lcmFreq
            indexflip = 1;
            %disp('over');
        end

        
    end 
      Priority(0); 
      frame_duration = Screen('GetFlipInterval', win);
    Screen('CloseAll');
    Screen('Close');
 
catch
    %sca
    Screen('CloseAll');
    Screen('Close');
    psychrethrow(psychlasterror);

end
end