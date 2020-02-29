% requirements: just Matlab

% 1.2.2020 Aljosa Klajderic, Hagebenerg

% brief: a Script that reads the legacy data fromt the previous project, trims the unrelevant values
% and plots.

clc; clear all; close all;

g = 9.81;
Fs = 400; % sampling rate
Ts = 1/Fs;
%%
% load data
load a_5g_26T.mat;
pjm   = a_5g_26T  * 1000 / g; % convert from m/s^2 to mg

pjm = resample(pjm,1,10); % resample pjm to 100 Hz


load LSM303_Crash3.mat; % LSM has freq of 400 Hz
xyz = LSM303_Crash3;
lms303 = xyz(:,2); % we are measuring only in y achse


%%

lms303 = lms303(1:end-100); % removing the old samples because something is wrong 

% trim the values that we do not need
[cbegin1, cend1]   = trimValues(lms303, 3, 150); % capture all the values bigger 110 mg and and additional capture 180 samples before and after
lms303    = lms303(cbegin1:cend1);



[cbegin2, cend2] = trimValues(pjm,2, 150);
pjm = pjm(cbegin2:cend2);

%SinglePlot2(lms303, 'Sensor LSM303', 1); % the one behind defines if you want to have another figure
%SinglePlot2(pjm, 'Sensor pjm', 0);

saveDataToTXT = 0;
if (saveDataToTXT == 1)
    
    fpathTXT = 'C:\Users\ak\Documents\Thesis_project\C++';
    
    fileID = fopen(fullfile(fpathTXT, 'lms303.txt'), 'w');
    fprintf(fileID,'%.4f\n', lms303);
    fclose(fileID);
    
    fileID = fopen(fullfile(fpathTXT, 'pjm.txt'), 'w');
    fprintf(fileID,'%.4f\n', pjm);
    fclose(fileID); 
end

size_pjm = size(pjm);
pjm_t = (1:size_pjm(1))*(1/Fs); % createa time vector for that

size_lms303 = size(lms303);
lms303_t = (1:size_lms303(1))*(1/Fs); % createa time vector for that


[mVal1, idx1] = maxk(lms303,2);



[mVal2, idx2] = maxk(pjm,2);
Integrate2 = trapz(pjm);
cdistance2 = cumtrapz(pjm);
size_cdistance2 = size(cdistance2);
cdistance_t2 = (0:size_cdistance2(1)-1)*(1/Fs); % createa time vector for that
T = table(cdistance_t2',pjm,cdistance2,'VariableNames',{'Time','Vektor','CumulativeDistance'});
disp(['first tbble with time: ', num2str(pjm_t(end)),' and max value: ',num2str(mVal2(1)),'\n']);

timepjm = 100;
plotBothSensorData(lms303,lms303_t - 0.0025,pjm, pjm_t, ' sensor lsm303', 'reference sensor pjm', timepjm, 0);
save=0;
if (save == 1 ) 
    fpath = 'C:\Users\ak\Desktop\\WIA3\git\IEEEtran\jpeg';
    saveas(gca, fullfile(fpath, 'ThirdCrash'), 'jpeg');
end % saving the fig at the specific location



%%
%functions
function plotBothSensorData(values1,time1,values2, time2, firstSens, SecSens, DurationOfThePulse, IntegrateVal)
    [mVal1, idx1] = maxk(values1,6);
    [mVal2, idx2] = maxk(values2,6);

    DoublePlot(values1,time1,values2, time2);
    xlim([0.0025 0.18]);
    %errorRate = abs(mVal1(2)-mVal2(2))/mVal1(2)*100;
    %str = ['Error rate for this signal is ', num2str(errorRate)];
    %title(['Time duration of the train impulse is ', num2str(DurationOfThePulse),' seconds, ',str, '%,','Integral value of Impulse is ',num2str(IntegrateVal) ,' mg']);
    title("Comparing acceleration data");
    hold on;
    
    plot(time1(idx1(1)),mVal1(1), 'y*','LineWidth',3 );
    plot(time1(idx1(2)),mVal1(2), 'y*','LineWidth',3 );
    plot(time1(idx1(3)),mVal1(3), 'y*','LineWidth',3 );
    plot(time1(idx1(4)),mVal1(4), 'y*','LineWidth',3 );
    plot(time1(idx1(5)),mVal1(5), 'y*','LineWidth',3 );
    %plot(time1(idx1(6)),mVal1(6), 'y*','LineWidth',3 );
    plot(time2(idx2(1)),mVal2(1), 'go', 'LineWidth',3);
    %plot(t_26T(idx2(2)),mVal2(2), 'go','LineWidth',3);
    %plot(t_26T(idx2(3)),mVal2(3), 'go','LineWidth',3);
    legend(firstSens,SecSens,...
            [' 1.Peak of 1. Signal at ', num2str(mVal1(1)),' mg'], ...
            [' 2.Peak of 1. Signal at ', num2str(mVal1(2)),' mg'], ...
            [' 3.Peak of 1. Signal at ', num2str(mVal1(3)),' mg'], ...
            [' 4.Peak of 1. Signal at ', num2str(mVal1(4)),' mg'], ...
            [' 5.Peak of 1. Signal at ', num2str(mVal1(5)),' mg'], ...
            [' 1.Peak of 2. Signal at ',  num2str(mVal2(1)),' mg']...
            );%'LOcation', 'southeast');
    
end

function pl = SinglePlot1(values, time,title1)
figure();
pl= plot(time,values);
hold on;
title(title1);
ylabel('acceleration [mg]');
xlabel('time [s]');
grid on;
%plot(values2);
end

function pl = SinglePlot2(values,title1, anotherFgure)
if anotherFgure == 1
    figure();
end;
pl = plot(values);%,'*');
hold on;
grid on;
title(title1);
legend(title1);
ylabel('acceleration [mg]');
xlabel('time [s]');
end

function DoublePlot1(values1,time1,values2, time2)
    figure();
    title('comparison between sensors');
    plot(time1, values1,'ro-'); % interpoliert
    hold on;
    grid on;
    plot(time2, values2, 'bh--');
    ylabel('acceleration [mg]'); % normal
    xlabel('time [s]');
    %hold off;
end

function DoublePlot(values1,time1,values2, time2)
    figure();
    title('comparison between sensors');
    plot(time1, values1,'r');
    hold on;
    grid on;
    plot(time2, values2, 'b');
    ylabel('acceleration [mg]');
    xlabel('time [s]');
    %hold off;
end
%{
    Function which trims values 
%}
function [cbegin,cend, cend1] = trimValues(vectorToTrim,margin,Limit)
    % saves values greater then Limit in idx and creates boolean array where condition is fullfiled 
    idx = vectorToTrim >=Limit; 
    
    % find gives us the indexes of the 1's 
    indexes = find(idx);
    
    % creating where we begin and end the vector
    cbegin = indexes(1)     - margin;
    cend   = indexes(end)   + margin;
    cend1  = indexes(end-1) + margin;
end



