
%   1)Commented out the 2 line block for importing data from .txt
%   2)Added input argument to function to be able to give existing array of BP, ECG,
%   3)Commented out the file selection 4 line block
%   4)Updated to remove reliance on 'filename' and 'pathname' variables from file selection when saving

function [final_data]=data_rr_detect_DLR(data_raw)

    
    % Perform QRS detection in order to obtain the RR intervals
    % Calibrate BP to unit mmHg
    % save the raw data in .mat file
    
    %{
    [filename, pathname] = uigetfile('*.txt','Select raw data file');                      % select the raw data file (Labview DAQ saves the data as .txt files)
    disp(['Processing ', filename]);
    %}
    
    prompt = {'Enter start time (in seconds)', 'Enter end time (in seconds)', 'Enter BP delay (in seconds)', 'Enter Height Correction (in cm)'};
    name = 'ENTER DATA SEGMENTATION INFORMATION';
    numlines = 1;
    answer = inputdlg(prompt, name, numlines); %--Dialog box

    
    fs = 1000;
    [b,a] = butter(4,20/(fs/2), 'low');           % Low Pass Filter for BP
    [b_l, a_l] = butter(4, 15/(fs/2), 'low');
    [b_h, a_h] = butter(4, 5/(fs/2), 'high');     % Band Pass Filter for ECG
    
    %{
    temp = importdata([pathname filename], '\t', 2);
    data_raw = temp.data;        % load raw data from .txt file
    %}
    if strcmp(answer{2},'end')
        answer{2} = num2str(floor(size(data_raw,1)/fs));
    end
    
    delay = round(str2double(answer{3})*fs);
    data = zeros(size(data_raw,1)-delay, size(data_raw,2));
    data(:,2) = data_raw(delay+1:end,2);
    data(:,[1 3:end]) = data_raw(1:end-delay,[1 3:end]);
    
    ind_s = str2double(answer{1})*fs+1;
    ind_e = str2double(answer{2})*fs;
    
    ecg = data(ind_s:ind_e,1);      % ECG: 1st channel
    bp = data(ind_s:ind_e,2).*100;  % BP: 2nd channel; 1V = 100mmHg (see manual of Finapres)
    data = data(ind_s:ind_e,:);
    
    hc = 1060*9.81*(str2double(answer{4})/100)*0.00750061683;
    bp = bp-hc;
    
    ecg = filtfilt(b_l, a_l, ecg);
    ecg = filtfilt(b_h, a_h, ecg);
    bp = filtfilt(b, a, bp);
    
    [qrs, qrsi, delay] = pan_tompkin_rev_xd(ecg,fs,0);    % obtain the time index for each QRS peak in the unit of sample (200Hz)
    %     qrsi = qrsi - delay;
    qrsi = qrsi.*(fs/200);       % convert to the sampling frequency (1000Hz)
    qrsi = qrsi(:);
    
    QRS = manudetect(qrsi, ecg, bp, fs);    % manually review
    
    % figure(1)
    % subplot(211); hold on; plot(bp); plot(ecg.*80+50,'r-'); stem(QRS, ecg(QRS).*80+50, 'g')
    % subplot(212); plot(QRS(2:end)./fs, diff(QRS).*1000./fs); xlabel('time (sec)'); ylabel('RR (ms)');
    
    [filenamesave, pathnamesave] = uiputfile('*.mat', 'Save as',[ '_' answer{1} '_' answer{2} '.mat']);
    
    save([pathnamesave filenamesave], 'bp', 'ecg', 'data', 'QRS', 'fs')
    
    final_data=load([pathnamesave filenamesave]); %create output structure
    disp('Done saving');
   
end