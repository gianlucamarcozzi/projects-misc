%%
clearvars

%% IMPORT
expName = "misc-E-002-" + ["006", "007", "010", "018"];
nExp = numel(expName);

% x-axis
mpfuAmp{1} = 0:0.1:99.9;
mpfuAmp{2} = 100:-0.1:0;
mpfuAmp{3} = 100:-0.2:0.1;
mpfuAmp{4} = 100:-0.2:0.1;
% Measuring I accidentally skipped some values (all the 0.9s for the first
% and all the 0.1s for the second)
idx1 = ones(size(mpfuAmp{1})) - [repmat([0, 1, zeros(1, 8)], 1, 100)];
mpfuAmp{1} = mpfuAmp{1}(logical(idx1));
idx2 = ones(size(mpfuAmp{2})) - [repmat([zeros(1, 8), 1, 0], 1, 100), 0];
mpfuAmp{2} = mpfuAmp{2}(logical(idx2));

for iexp = 1:nExp
    % Import data for channel 1 and channel 2
    if iexp < 3
        generalFolder = "../data/processed/";
        filename1 = generalFolder + expName(iexp) + "_ch1.txt";
        filename2 = generalFolder + expName(iexp) + "_ch2.txt";
    else
        generalFolder = "../data/raw/";
        filename1 = generalFolder + expName(iexp) + "-a.txt";
        filename2 = generalFolder + expName(iexp) + "-b.txt";
    end
    a1 = readtable(filename1);
    a2 = readtable(filename2);
    % Extract x and y values from the table
    y1{iexp} = table2array(a1(2:end, :));
    y2{iexp} = table2array(a2(2:end, :));
    x{iexp} = table2array(a1(1, :));
    % Number of measurements
    nMeas(iexp) = size(a2, 1) - 1;
    fprintf('Loaded %i out of %i\n', iexp, nExp)
end

%%
N_CH = 3;

% Define range of the x-axis where the pulses are expected
XRANGE{1} = [1.51, 1.695]*1e-6;
XRANGE{2} = [2.51, 2.695]*1e-6;
XRANGE{3} = [0.51, 0.695]*1e-6;
% Range of the baseline (outside)
XRANGEBL = [0.002, 3]*1e-6;

% Separate pulse 1 and pulse 2 from the y1 and y2 arrays
clear('ych', 'xch')
for iexp = 1:nExp
    for ich = 1:N_CH
        xCondition = x{iexp} > XRANGE{ich}(1) & x{iexp} < XRANGE{ich}(2);
        y1ch{iexp}{ich} = y1{iexp}(:, xCondition);
        y2ch{iexp}{ich} = y2{iexp}(:, xCondition);
        xch{iexp}{ich} = x{iexp}(xCondition);
    end
    xConditionBl = x{iexp} < XRANGEBL(1) | x{iexp} > XRANGEBL(2);
    bl1{iexp} = mean(y1{iexp}(:, xConditionBl), 2);
    bl2{iexp} = mean(y2{iexp}(:, xConditionBl), 2);
end

% Average
for iexp = 1:nExp
    % nm = nMeas(iexp);
    for ich = 1:N_CH
        y1fin{iexp}{ich} = mean(y1ch{iexp}{ich}, 2) - bl1{iexp};
        y2fin{iexp}{ich} = mean(y2ch{iexp}{ich}, 2) - bl2{iexp};
        amplitude{iexp}{ich} = hypot(y1fin{iexp}{ich}, y2fin{iexp}{ich});
        phase{iexp}{ich} = atan2(y2fin{iexp}{ich}, y1fin{iexp}{ich})*180/pi;
        
    end
end

 
%% --------------------------- PLOTS -------------------------------------
IEXP = 4;

figure(2)
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
for ii = 1:9
    nexttile
    xplot = x{IEXP};
    yplot = y1{IEXP}(ii, :);
    plot(xplot, yplot)
    hold on
    for ich = 1:N_CH 
        xplot = xch{IEXP}{ich};
        yplot = y1ch{IEXP}{ich}(ii, :);
        plot(xplot, yplot)
    end
    % ylim([0.008, 0.0095])
    % xlim(setaxlim(XRANGE{3}, 5))
end

figure(3)
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
for ii = 1:9
    nexttile
    xplot = x{IEXP};
    yplot = y2{IEXP}(ii, :);
    plot(xplot, yplot)
    hold on
    for ich = 1:N_CH 
        xplot = xch{IEXP}{ich};
        yplot = y2ch{IEXP}{ich}(ii, :);
        plot(xplot, yplot)
    end
    % ylim([0.008, 0.0095])
    % xlim(setaxlim(XRANGE{3}, 1.5))
end

% ------------------------------------------------------------------------
% -------------------------- FINAL FIGURE --------------------------------
% ------------------------------------------------------------------------
I2 = 2;
N_CH = 2;

figure(6)
clf
tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact")
nexttile(1)
hold on
box on
xplot = mpfuAmp{IEXP};
for ich = 1:N_CH
    yplot = y1fin{IEXP}{ich};
    plot(xplot, yplot, '.-')
end
xplot2 = mpfuAmp{I2};
for ich = 1:N_CH
    yplot2 = y1fin{I2}{ich};
    plot(xplot2, yplot2, 'Color', 'k')
end
title('channel 1')
xlim(setaxlim(xplot, 1))

nexttile(3)
hold on
box on
for ich = 1:N_CH
    yplot = y2fin{IEXP}{ich};
    plot(xplot, yplot, '.-')
end
for ich = 1:N_CH
    yplot2 = y2fin{I2}{ich};
    plot(xplot2, yplot2, 'Color', 'k')
end
title('channel 2')
xlim(setaxlim(xplot, 1))

nexttile(2)
hold on
box on
for ich = 1:N_CH
    plot(xplot, amplitude{IEXP}{ich}, '.-')
end
for ich = 1:N_CH
    plot(xplot2, amplitude{I2}{ich}, 'Color', 'k')
end
title('amplitude')
xlim(setaxlim(xplot, 1))

nexttile(4)
hold on
box on
for ich = 1:N_CH
    plot(xplot, phase{IEXP}{ich}, '.-')
end
for ich = 1:N_CH
    plot(xplot2, phase{I2}{ich}, 'Color', 'k')
end
title('phase')
xlim(setaxlim(xplot, 1))

% ------------------------------------------------------------------------
%------------------------- SAVE TO TXT -----------------------------------
% ------------------------------------------------------------------------
IEXP = 2;
chname = ["mpfu_BrX", "mpfu_BrMinX"];
for ich = 1:2
    savepath = chname(ich) + "_ampPhase_2025-03-28.txt";
    savemat = [mpfuAmp{IEXP}; amplitude{IEXP}{ich}'; phase{IEXP}{ich}'];
    % fid = fopen(savepath, 'w');
    % fprintf(fid, "%f %f %f\n", savemat);
    % fclose(fid);
end
