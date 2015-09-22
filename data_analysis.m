% rupert small, owlstone data analysis test
% 22.Sept.2015

clc
clear all
close all

data = csvread('test_matrix.csv',1);
peripheral = csvread('test_peripheral_dat.csv',1);
peripheral = peripheral(:,2:end); % remove row count column
[dm dn] = size(data);
[pm pn] = size(peripheral);
scanline = 1:pm;
norm_p = zeros(pm,pn);

% TASK 1 : investigate the peripheral data

% normalise each column
for i=1:pn
    norm_p(:,i) = (peripheral(:,i)-mean(peripheral(:,i)))/(max(peripheral(:,i))-min(peripheral(:,i)));
end
x_vals = 1:pn;
param_vals = 1:pm;
[X Y] = meshgrid(x_vals,param_vals);
figure

h=imagesc(norm_p);
set(gca, 'XTick', 1:pn)
set(gca, 'YTick', 2:2:pm)
title('peripheral data')
xlabel('feature')
ylabel('scan line number')
colorbar

% notable features
% 1. detector region temperature stays flat
figure
plot(scanline,peripheral(:,4),'b-o','MarkerFaceColor','r')
xlabel('scan line')
ylabel('detector region temp')
grid on

% 2. columns 7, 10, 12 and 13 are almost the same (after normalisation)
figure
hold on
plot(scanline,peripheral(:,7),'b-o','MarkerFaceColor','r')
plot(scanline,peripheral(:,10),'b-s','MarkerFaceColor','g')
plot(scanline,peripheral(:,12),'b-d','MarkerFaceColor','b')
plot(scanline,peripheral(:,13),'b-^','MarkerFaceColor','y')
xlabel('scan line')
ylabel('feature')
legend('inlet temp sht15', 'outlet temp','ctrl pcb temp','anlge pcb temp','Position',[200 70 200 70])
grid on

% 3. A smoothly varying column 6, with a flat partition
figure
plot(scanline,peripheral(:,6),'b-o','MarkerFaceColor','r')
xlabel('scan line')
ylabel('dispersion voltage')
grid on

% 4. The first column is binary - it takes one of two values
figure
[a b] = hist(peripheral(:,1),10);
bar(b, a)
xlabel('Time Step')
ylabel('dispersion voltage')
set(gca,'XTick',.001:.0001:.002)
grid on

% 5. Columns 2 and 14 are jagged / oscillatory
figure
subplot(2,1,1)
plot(scanline,peripheral(:,2),'b-o','MarkerFaceColor','r')
set(gca, 'XTickLabel','')
ylabel('diff pressure')
grid on
subplot(2,1,2)
plot(scanline,peripheral(:,14),'b-s','MarkerFaceColor','g')
xlabel('scan line')
ylabel('inlet flow')
grid on

% check if there's some periodicity there...(take FFT)
Y2 = fft(peripheral(:,2))/length(scanline);
Y14 = fft(peripheral(:,14))/length(scanline);
k = scanline(1:end-1)/max(scanline);
figure
subplot(2,1,1)
hold on
[ignore indices] = sort(abs(Y14(2:end)));
[ignore indices2] = sort(abs(Y2(2:end)));
for i=0:3
    x = k(indices(end-i));
    y = k(indices2(end-i));        
    rectangle('Position',[x,0,.02,.01],'FaceColor','y', 'LineStyle', 'none')
    rectangle('Position',[y,0,.02,.01],'FaceColor','y', 'LineStyle', 'none')
end
plot(k,abs(Y2(2:end)),'b-o','MarkerFaceColor','r')
set(gca, 'XTickLabel','')
ylabel('|FFT(diff pressure)|')
grid on
subplot(2,1,2)
hold on
for i=0:3
    x = k(indices(end-i));
    y = k(indices2(end-i));
    rectangle('Position',[x,0,.02,.04],'FaceColor','y', 'LineStyle', 'none')
    rectangle('Position',[y,0,.02,.04],'FaceColor','y', 'LineStyle', 'none')
end
plot(k,abs(Y14(2:end)),'b-s','MarkerFaceColor','g')
xlabel('freq^{-1}')
ylabel('|FFT(inlet flow)|')
grid on


% TASK 2 : investigate test_matrix
voltage = data(:,1);
faims = data(:,2:end);
figure
subplot(2,1,1)
imagesc(scanline,voltage,log(log(faims)))
set(gca,'YDir','normal')
set(gca, 'XTickLabel','')
ylabel('comp volt')
title('log-log FAIMS')
subplot(2,1,2)
contour(scanline,voltage,log(log(faims)))
xlabel('scanline')
ylabel('comp volt')
grid on

% find the FWHM for faims
voltage_delta = zeros(1, pm);
rownum = (1:dm)';
figure
hold on
for i=1:pm
    slice = data(:,i);
    size(slice)
    m = max(slice);
    if i<16
        [ignore indices] = sort(abs((slice - m/2)));
    else
        closest = min(abs(rownum - mm(1)), abs(rownum - mm(2)));
        [ignore indices] = sort(exp(closest.^6).*(abs(slice - m/2)));
    end
    subplot(18,3,i)
    plot(rownum(1:end-30),slice(1:end-30),'-k')
    set(gca, 'XTickLabel','','YTickLabel','','Box','off')
    xlim([0, max(rownum(1:end-30))])
    grid on
    mm = [min(indices(1:5)) max(indices(1:5))];
    voltage_delta(i) = voltage(mm(2)) - voltage(mm(1));
end

figure
bar(scanline, voltage_delta,'EdgeColor','b','FaceColor','w')
xlabel('scanline')
ylabel('$\Delta$ Voltage','Interpreter','LaTex')
xlim([0,max(scanline)])
grid on
