% requirements: just Matlab


% brief: a Script that connects to the oszi over TPC/IP and to the microcontroller
% over USB to received the acceleration from two different accelerometer on two different paths

%-----------------------------------------------------------
% Initialization:
%-----------------------------------------------------------
clc; clear all; close all;
instrfind;
instrreset;
disp("let's start :)");
convert_mg_s_ = 9.80665;
graph_move = 0;
to_add_to_offset = 1;
multipliction_constant = 5000;
try

    addpath('.\Oszi');
    rtb = VISA_Instrument('TCPIP::192.168.0.244::INSTR'); % Adjust the VISA Resource string to fit your instrument
    rtb.SetTimeoutMilliseconds(3000); % Timeout for VISA Read Operations
    
    idnResponse = rtb.QueryString('SYST:ERROR:ALL?'); % see if there are any errors in the error queue 
    rtb.ErrorChecking(); % Error Checking after Initialization block
    fprintf('\nvalue: %s\n', idnResponse);
    if ~(contains(idnResponse,"No error"))       
        exit(-1);
    end
       
    port = 'COM4';
    Baud_Rate = 115200;
    Data_Bits = 8;
    Stop_Bits = 1;
    %Default for Parity & Flow Control are none.  
    s = serial(port,'BaudRate',Baud_Rate,'DataBits',Data_Bits,'StopBits',Stop_Bits, 'Timeout', 25);
    s.Terminator = 'CR';
    fopen(s);
    disp("connected to port :)");
    
catch ME
        msgText = getReport(ME);
        disp(msgText);
        pause();
        fclose(s);
        rtb.Close() % Closing the session to the instrument 
end 

%-----------------------------------------------------------
% Main:
%-----------------------------------------------------------

try
    str = "y";
    while str ~= 'x' ||  isempty(str)
       rtb.Write('RUN');
       rtb.ErrorChecking(); % Error Checking after the data transfer
       disp('measurment started !') 
       while 1
         try
            idn = fscanf(s)
			% receiving the initial packet to get informed about the data size, threshold level, sampling time
            C = textscan(idn, '%s %u  %u %f %d \n');
            if contains(C{1}(1),"start")
               break;
            end
            pause(1);
         catch
            disp('timeout!')
         end  
       end    
       rtb.Write('STOP');
       rtb.ErrorChecking();


        HowLOngDidthePulseLast = C{3}(1); % cycles of the pulse
        dt = C{4}(1);  % sampling time
        threshold = double(C{5}) ;

        for time  = 1:(C{2})
            idn = fscanf(s);
            A = str2num(idn);
            x(time) = A(1) ;
            y(time) = A(2) ;
            z(time) = A(3) ;
            time;
        %time = time+1;
    
        end

        waveformBIN = rtb.QueryBinaryFloatData('FORM:BORD LSBF;:FORM REAL;:CHAN1:DATA?', false);
        samples_per_s = rtb.QueryDouble('ACQ:SRAT?'); % get the sample rate of the signal (samples/sec)
        rtb.ErrorChecking(); % Error Checking after the data transfer
        
        B = waveformBIN*multipliction_constant;
        
        % need to resample to 800 Hz exact; pjm sample time 4,17 k Hz 
        % (something is wrong oszi, could not change to any other lower sampling value that would make sens, only this was available)
        % 4,17k * 1000000 / 5208125 = 800
        if dt == 1/800
          B_new = resample(double(B), 1000000, 5208125);  % 800 Hz
        elseif dt == 1/400
          B_new = resample(double(B), 500000, 5208125);  % 400 Hz
        else
          B_new = resample(double(B), 250000, 5208125);  % 200 Hz
        end
        sample_time = (1/samples_per_s);
        size_oszi=size(B);
        size_oszi_new=size(B_new);
        t_oszi = (1:size_oszi(2))*(sample_time);
        t_oszi_new = (1:size_oszi_new(2))*(dt);
        

        
        z = -z;
        z = [zeros(1,20),z];
        x_size = size(z);
        t = (1:x_size(2))*dt;
        
        [max_oszi,idx_oszi] = max(B);
        [max_oszi_new,idx_oszi_new] = max(B_new);
        
        [max_lisdw12,idx_lisdw12] = max(z);
        
        graph_move= t_oszi(idx_oszi)-t(idx_lisdw12);
        graph_move_new= t_oszi_new(idx_oszi_new)-t(idx_lisdw12);
        
        end_lisdw12_time = t(end )- t(idx_lisdw12);
        
        
        bgn_cut_off_idx = find(abs(t_oszi_new-graph_move_new) < 0.001);
        end_cut_off_idx = find(abs(t_oszi_new-end_lisdw12_time) < 0.001);
     
        
        % cut off the vector to get the same length
        try
        B_new = B_new( (bgn_cut_off_idx(end)+to_add_to_offset): ( end_cut_off_idx(end) +idx_oszi_new ) );
        catch
            disp('something worng in the meas!!!')
            rtb.ErrorChecking(); 
            str = input("Continue\n",'s');
            continue;
        end
        % resample the values values
        size_oszi=size(B_new);
        t_oszi_new = (1:size_oszi(2))*(dt);
        %plot
        
        figure
        plot([0,t]+graph_move,[0,z], 'r-*');%, 'linewidth',2);
        hold on;
        plot([0,t_oszi],[0,B],'b--o'); % Displaying the waveform
        %plot(B); % Displaying the waveform

        ylabel('acceleration [mg]');
        legend("LISDW12","PJM LN");
        xlabel('time [s]');
        grid on;

        
        
        figure
        plot([0,t],[0,z], 'r-*');%, 'linewidth',2);
        hold on;
        plot([0,t_oszi_new],[0,B_new],'b--o'); % Displaying the waveform
        %plot(B); % Displaying the waveform

        ylabel('acceleration [mg]');
        legend("LISDW12","PJM LN");
        xlabel('time [s]');
        grid on;
        diff_new =  B_new(1:end)-z;
        plot(t,diff_new,'k:*'); % Displaying the waveform

            
        str = input("wating till the button is pressed\n",'s');
        if str == 'x'
            rtb.Write('RUN');
            break;
        end
        if str == 'y'
            if dt == 1/800 
              fpath = 'C:\Users\ak\Documents\Thesis_project\Matlab_Code\Data_Comparison\Dynamical_data_800Hz';
            elseif dt == 1/400 
              fpath = 'C:\Users\ak\Documents\Thesis_project\Matlab_Code\Data_Comparison\Dynamical_data_400Hz';
            else 
              fpath = 'C:\Users\ak\Documents\Thesis_project\Matlab_Code\Data_Comparison\Dynamical_data_200Hz';  
            end  
            Data.Lis2dw12.time = t;
            Data.Lis2dw12.timeoffset = graph_move;
            Data.Lis2dw12.x = x;
            Data.Lis2dw12.y = y;
            Data.Lis2dw12.z = z;
            
            Data.PJM.data = B_new;
            Data.PJM.time = t_oszi_new;
            
            Data.PJM.data_old = B;
            Data.PJM.time_old = t_oszi;
            
            Data.Lis2dw12.diff = diff_new;
            getTime = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
            save(fullfile(fpath, strcat('New-', getTime)),'Data')
            %saveas(gca, fullfile(fpath, getTime), 'jpeg');
        end 
        
        
        
        rtb.Write('RUN');
        rtb.ErrorChecking(); % Error Checking after the data transfer
        close all;
    end
    catch ME
        msgText = getReport(ME);
        disp(msgText);
         fclose(s);
         rtb.Close() % Closing the session to the instrument 
 
end
        fclose(s);
        rtb.ErrorChecking(); % Error Checking after the data transfer
        rtb.Close() % Closing the session to the instrument 
        %close all;
