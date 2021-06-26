% instal NI VISA 15.0 or newer
% This example does not require MATLAB Instrument Control Toolbox
% It is based on taking advantage of using .NET assembly called Ivi.Visa
% that is istalled together with NI VISA 15.0 or newer

% type "help VISA_Instrument" to get help on VISA_Instrument class


% brief: Script that in real times plots and refresh the data from the oszi and USB
clc; clear all; close all;

addpath('.\Oszi');
instrreset;
disp("let's start :)");

tx(1) = 0;
ty(1) = 0;
tz(1) = 0;
x(1) =  0;
y(1) = 0;
z(1) = 0;

mgForce(1) = 0;
t_oszi(1) = 0;

px = plot(tx,x, 'r');
hold on
py = plot(ty,y, 'b');%, 'linewidth',2);
pz = plot(tz,z, 'y');%, 'linewidth',2);
pB = plot(mgForce,t_oszi, ':');%, 'linewidth',2);
legend("x axis", "y axis","z axis");
title("Acceleration data");
ylabel('acceleration [mg]');
xlabel('time [s]');
%ylim([0 1000])
grid on;



px.YDataSource = 'x';
px.XDataSource = 'tx';
py.YDataSource = 'y';
py.XDataSource = 'ty';
pz.YDataSource = 'z';
pz.XDataSource = 'tz';
pz.YDataSource = 'z';
pz.XDataSource = 'tz';
pB.YDataSource = 'mgForce';
pB.XDataSource = 't_oszi';
MeanValue = 0.03279;
try
    port = 'COM4';
    Baud_Rate = 115200;
    Data_Bits = 8;
    Stop_Bits = 1;
    %Default for Parity & Flow Control are none.  
    s = serial(port,'BaudRate',Baud_Rate,'DataBits',Data_Bits,'StopBits',Stop_Bits, 'Timeout', 12);
    s.Terminator = 'CR';
    fopen(s);
    disp("connected to port :)");
    
    
    rtb = VISA_Instrument('TCPIP::192.168.0.244::INSTR'); % Adjust the VISA Resource string to fit your instrument
    rtb.SetTimeoutMilliseconds(3000); % Timeout for VISA Read Operations
    
    idnResponse = rtb.QueryString('SYST:ERROR:ALL?'); % see if there are any errors in the error queue 
    rtb.ErrorChecking(); % Error Checking after Initialization block
    fprintf('connected to Oszi error: %s\n', idnResponse);
    
    pause on;

    disp('press the button :D');  
    time = -1;
    Mean = zeros(30, 2,'double');
    Meanidx = 1;
    while Meanidx<=30

    idn = fscanf(s);
    C = textscan(idn, '%s %u  %u %f %d \n');
    dt = C{4}(1);  % sampling time
    num_of_elems = C{2};
    tx = zeros(1, num_of_elems,'double');
    ty = zeros(1, num_of_elems,'double');
    tz = zeros(1, num_of_elems,'double');
    x = zeros(1, num_of_elems,'double');
    y = zeros(1, num_of_elems,'double');
    z = zeros(1, num_of_elems,'double');
    y = zeros(1, num_of_elems,'double');
    z = zeros(1, num_of_elems,'double');
    
    time = 1;
    idn={''};
    
    for time  = 1:(C{2})
        
        idn = fscanf(s);
        A = str2num(idn);
        x(time) = A(1)  ;
        y(time) = A(2)  ;
        z(time) = A(3)*(-1);
        tx(time) = double(time)*dt;
        %tx(time)
        ty(time) = tx(time);
        tz(time) = ty(time);
        time = time+1;
    end
  
    refreshdata;
    drawnow;

    Meanidx = Meanidx + 1;
   end
catch ME
   fclose(s);
   rtb.Close(); % Closing the session to the instrument 
   msgText = getReport(ME);
   disp(msgText);
end
Mean60_2 = Mean;
save('_60Degress_2.mat','Mean60_2')
 fclose(s);
 rtb.Close() % Closing the session to the instrument 
 
