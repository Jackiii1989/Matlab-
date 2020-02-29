% instal NI VISA 15.0 or newer
% This example does not require MATLAB Instrument Control Toolbox
% It is based on taking advantage of using .NET assembly called Ivi.Visa
% that is istalled together with NI VISA 15.0 or newer

% type "help VISA_Instrument" to get help on VISA_Instrument class

% 1.2.2020 Aljosa Klajderic, Hagebenerg

% brief: a Script that connects to the oszi over TPC  witht VISA driver and reads 
% the current measurment with the configured time base in oszi

%-----------------------------------------------------------
% Initialization:
%-----------------------------------------------------------
instrreset;
clc;clear all; close all
try
    rtb = VISA_Instrument('TCPIP::192.168.0.244::INSTR'); % Adjust the VISA Resource string to fit your instrument
    rtb.SetTimeoutMilliseconds(3000); % Timeout for VISA Read Operations
catch ME
    error ('Error initializing the instrument:\n%s', ME.message);
end

try
    idnResponse = rtb.QueryString('*IDN?'); % query the instruments identificstion string 
    fprintf('\nInstrument Identification string: %s\n', idnResponse);
   
   % Reset the instrument, clear the Error queue and sanity check if there are any errors 
    rtb.Write('*RST;*CLS'); 
	idnResponse = rtb.QueryString('SYST:ERROR:ALL?'); % see if there are any errors in the error queue 
	fprintf('\nError value: %s\n', idnResponse);

%-----------------------------------------------------------
% Main:
%-----------------------------------------------------------

    str = "y";
    while str ~= 'x' ||  isempty(str)
        close all;
        
        fprintf('Fetching waveform in binary format... ');
        waveformBIN = rtb.QueryBinaryFloatData('FORM:BORD LSBF;:FORM REAL;:CHAN1:DATA?', false);
        fprintf('Samples count: %d\n', size(waveformBIN, 2));
        size_oszi=size(waveformBIN);
        samples_per_s = rtb.QueryDouble('ACQ:SRAT?') % get the sample rate of the signal (samples/sec)
        channelOffset = rtb.QueryDouble('CHAN1:OFFS?');
        rtb.ErrorChecking(); % Error Checking after the data transfer
        sample_time = (1/samples_per_s);%/(25*10e10);
        t_oszi = (1:size_oszi(2))*(sample_time);
        mean(waveformBIN)
        aB = waveformBIN - mean(waveformBIN);
        plot(t_oszi,aB); % Displaying the waveform
        grid on;
        ylim([-0.1 0.1]);
        pause(0.1);
        str = input("wating till the button is pressed\n",'s');
        %str = 'x';
        
    end
    
    rtb.Close() % Closing the session to the instrument 
    % look if they are any static errors after closing
    %idnResponse = rtb.QueryString('SYST:SERR?'); % see if there are any errors in the error queue 
    %rtb.ErrorChecking(); % Error Checking after Initialization block
    %fprintf('\nError value: %s\n', idnResponse);
    
    % -----------------------------------------------------------
    % Error handling
    % -----------------------------------------------------------
catch ME
    switch ME.identifier
        case 'VISA_Instrument:ErrorChecking'
            rethrow(ME);
        otherwise
            rethrow(ME);
    end
end