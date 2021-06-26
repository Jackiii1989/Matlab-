%clc; 
clear all; close all;
                                             
% define variables 
dinfo = dir('new-2019*.mat');
fpath = 'C:\Documents\Thesis_project\Matlab_Code\Data_Comparison';
save = 0;
    
Begin_iter = 1;
End_iter = length(dinfo);

figuredisp = [0 0 0 1 1];

x_y_vector = zeros( End_iter, 4);
if figuredisp(3) == 1
    f_dif = figure('units','normalized','outerposition',[0 0 1 1]);
else
    f_dif = 0;
end

for i = Begin_iter:End_iter
   
    load(dinfo(i).name);
    [ pjm, diff, plot_numb, Integrate_error] = analyze( Data, f_dif, figuredisp, i, fpath, save);
    x_y_vector(i,1)  = pjm;
    x_y_vector(i,2)  = diff;
    x_y_vector(i,3)  = i; % the index of the dinfo to know to which the files belongs to. 
    x_y_vector(i,4)  = plot_numb; % the index of the dinfo to know to which the files belongs to. 
    x_y_vector(i,5)  = Integrate_error;  
    
end

if figuredisp(3) == 1
    figure(f_dif);

    fpath = 'C:\Users\aljosa.klajderic\Documents\Thesis_project\Matlab_Code\Data_Comparison';
    
    if save == 1
    saveas(gca, fullfile(fpath, 'diff_data_800'), 'png');
    end % saving the fig at the specific location
end

accumalate_peak = 0;
accumalate_integ = 0;
x_y_vector_sorted = sortrows(x_y_vector);
for x = 1:length(x_y_vector_sorted(:,2))
    accumalate_peak = accumalate_peak + abs(x_y_vector_sorted(x,2));
    accumalate_integ = accumalate_integ + abs(x_y_vector(x,5));
end

res_peak = accumalate_peak / length(x_y_vector_sorted(:,2));
res_integ = accumalate_integ / length(x_y_vector(:,5));

if figuredisp(4) == 1
    figure('units','normalized','outerposition',[0 0 1 1])

    title("Comparison data at 800 khz");

    plot(x_y_vector_sorted(:,1),abs(x_y_vector(:,2)),'--bo'); % Displaying the waveform
    hline = refline([0 res_peak]);
    hline.Color = 'r';

    legend("lisdw12 error peak rate", "average peak error rate");
    ylabel('Absolute error values of lisd12[mg]');
    xlabel('Amplitude values[mg]');
    xlim([x_y_vector_sorted(1,1) x_y_vector_sorted(end,1)] );
    ylim([-10 1300] );
    grid minor;
    grid on;

    for i= 1:length(x_y_vector_sorted(:,1))
        text(x_y_vector_sorted(i,1), abs(x_y_vector(i,2)),int2str(i), 'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
end

set(gca,'FontSize',20)
if save == 1 

    saveas(gca, fullfile(fpath, 'error_peak_800'), 'png');
end % saving the fig at the specific location

%figure
    if  figuredisp(5) == 1
    figure('units','normalized','outerposition',[0 0 1 1])
    % error integration rate
    plot(x_y_vector_sorted(:,1),abs(x_y_vector(:,5)),'--bo'); % Displaying the waveform
    hline = refline([0 res_integ]);
    hline.Color = 'r';

    for i= 1:length(x_y_vector_sorted(:,1))
        text(x_y_vector_sorted(i,1), abs(x_y_vector(i,5)),int2str(i), 'VerticalAlignment','bottom','HorizontalAlignment','right')
    end

    legend("error integration rate of the lisdw12", "average integration error rate amongst all measurements");
    ylabel('Absolute error values of lisd12[mg]');
    xlabel('Amplitude values[mg]');
    xlim([x_y_vector_sorted(1,1) x_y_vector_sorted(end,1)] );
    ylim([0 7] );
    grid minor;
    grid on;
    set(gca,'FontSize',20)

    if save == 1 
        saveas(gca, fullfile(fpath, 'error_integ_800'), 'png');
    end % saving the fig at the specific location
end
disp(['average error value of the peak is ', num2str(res_peak)]);
disp(['average error value of the peak is ', num2str(res_integ)]);
disp('code finshed for 800 Hz!!!')


function [ max_pjm, diff, plot_numb, Integrate_error] = analyze( Data, f_dif, figuredisp, i, fpath, save)

    fft_enable = 0;
    cut_off_enable= 0;
    plot_numb = 0 ;
    convert_mg_s_ = 9.80665/1000;
    lsdw12_t   = Data.Lis2dw12.time;
    lsdw12_z   = Data.Lis2dw12.z;
    graph_move = Data.Lis2dw12.timeoffset;

    pjm_resample = Data.PJM.data;
    pjm_t_resample = Data.PJM.time;
    
    pjm   = Data.PJM.data_old;
    pjm_t = Data.PJM.time_old;
    diff_vec = Data.Lis2dw12.diff;
        
    if figuredisp(1) == 1
        %figure
        figure('units','normalized','outerposition',[0 0 1 1])
        %set(gcf,'position',[10 490 603 500])
        
        
        hold on;
        diff47=lsdw12_t(2)-lsdw12_t(1);
        plot(lsdw12_t + graph_move,lsdw12_z-diff47, 'r-*');%, 'linewidth',2);
        plot(pjm_t,pjm,'b--o'); % Displaying the waveform
        legend("lsd2dw12 data", "pjm  data");
        title(['Comparing raw data of measurement ', int2str(i)]);
        ylabel('acceleration [mg]');
        xlabel('time [s]');
        xlim([graph_move+lsdw12_t(10) lsdw12_t(30)+graph_move] );
        ylim([-2000 5500] );
        set(gca,'Box','on');
        grid on;
        hold off;
        set(gca,'FontSize',20)
        if save == 1
            %figure('units','normalized','outerposition',[0 0 1 1])
            name = strcat('raw-data-800Hz', int2str(i));
            saveas(gca, fullfile(fpath, name), 'png');
        end % saving the fig at the specific location

        
    end    
    
    
    
    
    if figuredisp(2) == 1
        %figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        %[left bottom width height]
        %set(gcf,'position',[1230 490 603 500])
        %set(gcf,'position',[620 280 603 500])
        figure('units','normalized','outerposition',[0 0 1 1])
        plot(lsdw12_t,lsdw12_z, 'r-*');%, 'linewidth',2);
        hold on;
        plot(pjm_t_resample,pjm_resample,'b--o'); % Displaying the waveform
        
        %plot(pjm_t_resample,diff_vec,'k:d'); % Displaying the waveform
        title(['Comparing resampled data of measurement ', int2str(i)]);
        legend("lsd2dw12 data", "pjm resampled data");
        %legend("lsd2dw12 data", "pjm resampled data", "amplitude difference between lsd2dw12 and pjm");
        ylabel('acceleration [mg]');
        xlabel('time [s]');
        ylim([-2000 4000] );
        xlim([lsdw12_t(18) lsdw12_t(30)] );
        %title(name);
        grid on;
        set(gca,'FontSize',20);
        if cut_off_enable ==1
            xlim(time_cut_off{1} );  
        end
        
        if fft_enable == 1
          figure
          L = size(pjm_t_resample);
          Y = fft(pjm_resample);
          P1 = abs(Y/L(2));
          Fs = 1/pjm_t_resample(2)-pjm_t_resample(1);
         %f = Fs*(0:(L(2)/2))/L(2);
          f = Fs*(0:(L(2)))/L(2);
          plot(f(1:end-1),P1) 
          title('Single-Sided Amplitude Spectrum of X(t)')
          xlabel('f (Hz)')
          ylabel('|P1(f)|')
          grid on;
          
        end
        if save == 1
            name = strcat('resampled-data-at-800-Hz', int2str(i));
            saveas(gca, fullfile(fpath, name), 'png');
        end % saving the fig at the specific location
        
    end
    
    % return the x y coordinated
    [max_lsdw12, idx_lsdw12] = max(lsdw12_z);
    [max_pjm, idx_new_pjm] = max(pjm_resample);
    diff = max_pjm  - max_lsdw12;
    
    %graph_move= new_pjm_t(idx_new_pjm)-lsdw12_t(idx_lsdw12);
    
    
    
    if figuredisp(3) == 1
        figure(f_dif);
         %[left bottom width height]
        %set(gcf,'position',[2450 490 603 500])
        hold on;
        plot_numb = plot(pjm_t_resample,diff_vec); % Displaying the waveform
        title("Differnce Values between two accelerometers");
        ylabel('acceleration [mg]');
        xlabel('time [s]');
        ylim([-2000 3000]);
        xlim([0 0.12]);
        grid on;     
    end

    % find the first negative value in the impact
    idx_end_integrate_lsdw12 = find(lsdw12_z<0, 1);
    
    
    idx_begin_integrate_pjm = 0;
    for x = 0:20
        if pjm_resample(idx_new_pjm-x)<0
            idx_begin_integrate_pjm = idx_new_pjm-x;
            break;
        end
    end
    
    idx_end_integrate_pjm  = 0;
    for x = idx_new_pjm:length(pjm_resample)
        if (pjm_resample(x)< 0.01)
            idx_end_integrate_pjm = x;
            break;
        end
    end

    
    
    lsdw12_z_integrate = lsdw12_z(20:idx_end_integrate_lsdw12);
    pjm_resample_integrate = pjm_resample((idx_begin_integrate_pjm):idx_end_integrate_pjm);
    
    integration_lsdw12_t = pjm_t_resample(20:idx_end_integrate_lsdw12);
    integration_pjm_t = pjm_t_resample((idx_begin_integrate_pjm):idx_end_integrate_pjm);
    %length(integration_c)
    
    Integrate_lisdw12 = trapz(integration_lsdw12_t, lsdw12_z_integrate);%*convert_mg_s_;
    Integrate_pjm     = trapz(integration_pjm_t, pjm_resample_integrate);%*convert_mg_s_;
    Integrate_error =  Integrate_pjm - Integrate_lisdw12;
    %diff_vector = pjm-lsdw12_z; 
    
    
end
