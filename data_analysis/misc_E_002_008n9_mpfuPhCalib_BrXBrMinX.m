%%
clearvars

%% IMPORT
generalFolder = "../data/raw/";
expName = ["misc-E-002-008", "misc-E-002-009"];
nExp = numel(expName);

for iexp = 1:nExp
    % Import phases from the param file
    apar = readtable(generalFolder + expName(iexp) + "-param.txt", ...
        'ImportErrorRule', 'omitrow');
    mpfuPh{iexp} = table2array(apar(:, 2));
    % Import data for channel 1 and channel 2
    a1 = readtable(generalFolder + expName(iexp) + ".txt");
    % Extract x and y values from the table
    y1{iexp} = table2array(a1(2:end, :));
    x{iexp} = table2array(a1(1, :));
    % Number of measurements
    nMeas(iexp) = size(a1, 1) - 1;
end

%%
N_CH = 3;
IEXP = 2;

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
        xch{iexp}{ich} = x{iexp}(xCondition);
    end
    xConditionBl = x{iexp} < XRANGEBL(1) | x{iexp} > XRANGEBL(2);
    bl1{iexp} = mean(y1{iexp}(:, xConditionBl), 2);
end

% Average
for iexp = 1:nExp
    for ich = 1:N_CH
        y1fin{iexp}{ich} = mean(y1ch{iexp}{ich}, 2) - bl1{iexp};
    end
end

%% ---------------------------- PLOTS ------------------------------------
% ------------------------------------------------------------------------
figure(2)
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
for ii = 1:9
    nexttile
    hold on
    box on
    for ich = 1:N_CH 
        xplot = xch{IEXP}{ich};
        yplot = y1ch{IEXP}{ich}(ii, :);
        plot(xplot, yplot)
    end
    xplot = x{IEXP};
    yplot = y1{IEXP}(ii, :);
    h0 = plot(xplot, yplot, 'k');
    uistack(h0, 'bottom')
    % ylim([0.008, 0.0095])
    % xlim(setaxlim(xplot, 1))
end

% ------------------------------------------------------------------------
% -------------------------- FINAL FIGURE --------------------------------
% ------------------------------------------------------------------------
figure(4)
clf

xxs = [17.5, 43, 74];
xx1 = linspace(xxs(1), xxs(2), 100);
xx2 = linspace(xxs(2), xxs(3), 100);

tL = tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");
nexttile
hold on
box on
xplot = mpfuPh{IEXP};
for ich = 1:2
    plot(xplot, y1fin{IEXP}{ich}, 'o-')
    % plot(xplot + (ich - 2)*1, y1fin{IEXP}{ich}, 'o-')
end
xline(xxs)
yline(0, 'Color', [0.8, 0.8, 0.8], 'HandleVisibility', 'off')
xlim(setaxlim(xplot, 1))
yticks(-0.6:0.3:0.6)
legend('<+x>', ['<', char(hex2dec('2013')), 'x>'], 'Location', 'north')
%
nexttile
hold on
box on
for iexp = 1:nExp
    pulsePh = acos(y1fin{iexp}{ich}/max(abs(y1fin{iexp}{ich})));
    plot(xplot, pulsePh, '.-')
end
hold on
plot(xx1, xx1*7.*(pi/180) - 2.15, 'r')
plot(xx2, xx2*(-5.7)*(pi/180) + 7.3, 'r')
xline(xxs)
labelaxesfig(tL, 'MPFU phase / %', "TM Phase / a.u.")

%% 
figure(1)
clf
hold on
box on

IEXP = 1;
xplot = mpfuPh{IEXP};
for ich = 1:2
    plot(xplot + (ich-1) * 1.44, ...
        y1fin{IEXP}{ich}/max(abs(y1fin{IEXP}{ich})), '.-')
    % plot(xplot + (ich - 2)*1, y1fin{IEXP}{ich}, 'o-')
end

IEXP = 2;
xplot = mpfuPh{IEXP};
for ich = 1:2
    plot(xplot + (ich-1) * 0.5 + 1.5, y1fin{IEXP}{ich}/max(abs(y1fin{IEXP}{ich})), '.-')
    % plot(xplot + (ich - 2)*1, y1fin{IEXP}{ich}, 'o-')
end

IEXP = 1;
actualPhase = acos(y1fin{IEXP}{ich}/max(abs(y1fin{IEXP}{ich})));
figure(10)
plot(xplot, actualPhase)
hold on
plot(xplot, ...
        y1fin{IEXP}{ich}/max(abs(y1fin{IEXP}{ich})), '.-')

%% ------------------------- SAVE TO TXT ---------------------------------
IEXP = 1;
chname = ["MpfuBrX", "MpfuBrMinX"];
for ich = 1:2
    savepath = "calib" + chname(ich) + "_2025-03-28.txt";
    savemat = [mpfuPh{IEXP}; amplitude{IEXP}{ich}'; phase{IEXP}{ich}'];
    % fid = fopen(savepath, 'w');
    % fprintf(fid, "%f %f %f\n", savemat);
    % fclose(fid);
end
