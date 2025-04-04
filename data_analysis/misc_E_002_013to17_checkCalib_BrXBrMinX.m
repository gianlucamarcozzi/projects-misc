%%
clearvars

%% IMPORT
generalFolder = "../data/raw/";
% generalFolder = "/home/gianluca/files/mnt/1/files/projects/misc/misc-E-002-012/";
expName = ["misc-E-002-013", "misc-E-002-014", "misc-E-002-015", ...
    "misc-E-002-016", "misc-E-002-017"];
nExp = numel(expName);

% x-axis
pulseAmps{1} = linspace(0.06, 0.0014, 10)./0.06;
pulseAmps{2} = linspace(0.06, 0.0014, 25)./0.06;
pulseAmps{3} = linspace(0.06, 0.0014, 25)./0.06;
pulseAmps{4} = linspace(0.06, 0.0014, 25)./0.06;
pulseAmps{5} = linspace(0.06, 0.0014, 25)./0.06;

for iexp = 1:nExp
    % Import data for channel 1 and channel 2
    a1 = readtable(generalFolder + expName(iexp) + "-a.txt");
    a2 = readtable(generalFolder + expName(iexp) + "-b.txt");
    % Extract x and y values from the table
    y1{iexp} = table2array(a1(2:end, :));
    y2{iexp} = table2array(a2(2:end, :));
    x{iexp} = table2array(a1(1, :));
    % Number of measurements
    nMeas(iexp) = size(a2, 1) - 1;
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
    for ich = 1:N_CH
        y1fin{iexp}{ich} = mean(y1ch{iexp}{ich}, 2) - bl1{iexp};
        y2fin{iexp}{ich} = mean(y2ch{iexp}{ich}, 2) - bl2{iexp};
        amplitude{iexp}{ich} = hypot(y1fin{iexp}{ich}, y2fin{iexp}{ich});
        phase{iexp}{ich} = atan2(y2fin{iexp}{ich}, y1fin{iexp}{ich})*180/pi;
    end
end

%% --------------------------- PLOTS -------------------------------------
% ------------------------------------------------------------------------
IEXP = 5;

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
    % xlim(setaxlim(XRANGE{2}, 5))
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
    % xlim(setaxlim(XRANGE{2}, 1.5))
end

% ------------------------------------------------------------------------
% -------------------------- FINAL FIGURE --------------------------------
% ------------------------------------------------------------------------
% IEXP2 = 4;
N_CH = 2;
figure(55)
clf
tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact")
nexttile(1)
hold on
box on
xplot = pulseAmps{IEXP};
for ich = 1:N_CH
    plot(xplot, y1fin{IEXP}{ich}, 'o-')
end
% xplot2 = pulseAmps{IEXP2};
% for ich = 1:N_CH
%     plot(xplot2, y1fin{IEXP2}{ich}, 'o-')
% end
title('channel 1')
xlim(setaxlim(xplot, 1))

nexttile(3)
hold on
box on
for ich = 1:N_CH
    plot(xplot, y2fin{IEXP}{ich}, 'o-')
end
% for ich = 1:N_CH
%     plot(xplot2, y2fin{IEXP2}{ich}, 'o-')
% end
title('channel 2')
xlim(setaxlim(xplot, 1))

nexttile(2)
hold on
box on
for ich = 1:N_CH
    plot(xplot, amplitude{IEXP}{ich}, 'o-')
end
% for ich = 1:N_CH
%     plot(xplot2, amplitude{IEXP2}{ich}, 'o-')
% end
title('amplitude')
xlim(setaxlim(xplot, 1))

nexttile(4)
hold on
box on
for ich = 1:N_CH
    plot(xplot, 180*(ich-1) + phase{IEXP}{ich}, 'o-')
end
% for ich = 1:N_CH
%     plot(xplot2, 180*(ich-1) + phase{IEXP2}{ich}, 'o-')
% end
title('phase')
xlim(setaxlim(xplot, 1))

% ------------------------------------------------------------------------
%------------------------- SAVE TO TXT -----------------------------------
% ------------------------------------------------------------------------
IEXP = 2;
chname = ["mpfu_BrX", "mpfu_BrMinX"];
for ich = 1:2
    savepath = chname(ich) + "_ampPhase_2025-03-28.txt";
    savemat = [mpfuAmps; amplitude{IEXP}{ich}'; phase{IEXP}{ich}'];
    % fid = fopen(savepath, 'w');
    % fprintf(fid, "%f %f %f\n", savemat);
    % fclose(fid);
end
