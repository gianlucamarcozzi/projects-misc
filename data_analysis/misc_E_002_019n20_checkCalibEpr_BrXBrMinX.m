%%
clearvars

%% IMPORT
generalFolder = "../data/raw/";
% generalFolder = "/home/gianluca/files/mnt/1/files/projects/misc/misc-E-002-012/";
expName = ["misc-E-002-020"];
nExp = numel(expName);

% x-axis
pulseAmps = linspace(0.06, 0.0014, 25)./0.06;

for iexp = 1:nExp
    expFolder = generalFolder + expName(iexp);
    % BrX and BrMinX standing echo data (p plus, m minus)
    [xp{iexp}, yp{iexp}, parp{iexp}] = ...
        loadfolderelexsys(expFolder, '/*a.DTA');
    [xm{iexp}, ym{iexp}, parm{iexp}] = ...
        loadfolderelexsys(expFolder, '/*b.DTA');
    % Import data for channel 1 and channel 2 of the oscilloscope
    a1 = readtable(expFolder + "/" + expName(iexp) + "-scope-ch1.txt");
    a2 = readtable(expFolder + "/" + expName(iexp) + "-scope-ch2.txt");
    % Extract x and y values from the table
    ych1{iexp} = table2array(a1([2, 4:end], :));
    ych2{iexp} = table2array(a2([2, 4:end], :));
    xscope{iexp} = table2array(a1(1, :));
    % Number of measurements
    nMeas(iexp) = (size(a2, 1) - 2)/2;
end
xx{iexp} = xp{iexp}{1};

isc = 0;
for iexp = 1:nExp
    for ii = 1:2:size(ych1{iexp}, 1)
        isc = isc + 1;
        ysc1{iexp}(isc, :) = ych1{iexp}(ii,:) + 1i*ych2{iexp}(ii,:);
        ysc2{iexp}(isc, :) = ych1{iexp}(ii+1,:) + 1i*ych2{iexp}(ii+1,:);
    end
end

%% OSCILLOSCOPE
N_PULSE = 2;

% Define range of the x-axis where the pulses are expected
XRANGE{1} = [3.05, 3.10]*1e-7;
XRANGE{2} = [7.06, 7.20]*1e-7;
% Range of the baseline (outside)
XRANGEBL = 10*1e-7;

% Separate pulse 1 and pulse 2 from the y1 and y2 arrays
clear('ych', 'xch')
for iexp = 1:nExp
    for ip = 1:N_PULSE
        xCondition = xscope{iexp} > XRANGE{ip}(1) & ...
            xscope{iexp} < XRANGE{ip}(2);
        y1ch{iexp}{ip} = ysc1{iexp}(:, xCondition);
        y2ch{iexp}{ip} = ysc2{iexp}(:, xCondition);
        xch{iexp}{ip} = xscope{iexp}(xCondition);
    end
    xConditionBl = xscope{iexp} > XRANGEBL;
    bl1{iexp} = mean(ysc1{iexp}(:, xConditionBl), 2);
    bl2{iexp} = mean(ysc2{iexp}(:, xConditionBl), 2);
end

% Average
for iexp = 1:nExp
    for ip = 1:N_PULSE
        y1fin{iexp}{ip} = mean(y1ch{iexp}{ip}, 2) - bl1{iexp};
        y2fin{iexp}{ip} = mean(y2ch{iexp}{ip}, 2) - bl2{iexp};
        amplitude1{iexp}{ip} = abs(y1fin{iexp}{ip});
        amplitude2{iexp}{ip} = abs(y2fin{iexp}{ip});
        phase1{iexp}{ip} = ... 
            atan2(imag(y1fin{iexp}{ip}), real(y1fin{iexp}{ip}))*180/pi;
        phase2{iexp}{ip} = ... 
            atan2(imag(y2fin{iexp}{ip}), real(y2fin{iexp}{ip}))*180/pi;
    end
end
%% ------------------------ PLOTS SCOPE ----------------------------------
% ------------------------------------------------------------------------
IEXP = 1;

figNames = ["BrX real", "BrX imag", "BrMinX real", "BrMinX imag"];
funplots = split(figNames(:), ' ');
funplots = funplots(:, 2);  % real, imag, real, imag
for iplot = 1:4
    f = figure(iplot + 1);  % 2 to 4
    f.Name = figNames(iplot);
    funplot = funplots(iplot);
    clf
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
    for ii = 1:9
        nexttile
        xplot = xscope{IEXP};
        if iplot < 3
            yplot = ysc1{IEXP}(ii, :);
        elseif iplot < 5
            yplot = ysc2{IEXP}(ii, :);
        else
            error('Error with plot.')
        end
        yplot = feval(funplot, yplot);
        plot(xplot, yplot)
        hold on
        for ip = 1:N_PULSE 
            if iplot < 3
                yplot = y1ch{IEXP}{ip}(ii, :);
            elseif iplot < 5
                yplot = y2ch{IEXP}{ip}(ii, :);
            else
                error('Error with plot.')
            end
            xplot = xch{IEXP}{ip};
            yplot = feval(funplot, yplot);
            plot(xplot, yplot)
        end
    end
end

%% ------------------- FINAL FIGURE SCOPE --------------------------------
figure(11)
axTitles = ["TM real", "TM imag"];
funplots = split(axTitles(:), ' ');
funplots = funplots(:, 2);  % real, imag
clf
tiledlayout(2, 3, "TileSpacing", "compact", "Padding", "compact")
% Plots 1 and 3: real and imag part
numTile = [1, 4];
for iplot = 1:2
    nexttile(numTile(iplot))
    hold on
    box on

    xplot = pulseAmps;
    plot(xplot, feval(funplots(iplot), y1fin{IEXP}{1}), 'o-')  % <+x>
    plot(xplot, feval(funplots(iplot), y2fin{IEXP}{1}), 'o-')  % <-x>
    plot(xplot, feval(funplots(iplot), y1fin{IEXP}{2}), 'o-')  % +x pi
    plot(xplot, feval(funplots(iplot), y2fin{IEXP}{2}), 'o-')  % +x pi

    title(axTitles(iplot))
    xlim(setaxlim(xplot, 1))
end

% Amplitude
nexttile(2)
hold on
box on
plot(xplot, amplitude1{IEXP}{1}, 'o-')  % <+x>
plot(xplot, amplitude2{IEXP}{1}, 'o-')  % <-x>
plot(xplot, amplitude1{IEXP}{2}, 'o-')  % +x pi
plot(xplot, amplitude2{IEXP}{2}, 'o-')  % +x pi (same)
title('Amplitude')
xlim(setaxlim(xplot, 1))

% Phase
nexttile(3)
hold on
box on
plot(xplot, phase1{IEXP}{1}, 'o-')  % <+x>
plot(xplot, phase2{IEXP}{1}, 'o-')  % <-x>
plot(xplot, phase1{IEXP}{2}, 'o-')  % +x pi
plot(xplot, phase2{IEXP}{2}, 'o-')  % +x pi (same)
title('phase')
xlim(setaxlim(xplot, 1))

% Relative amplitude difference
nexttile(5)
hold on
box on
yplot = (amplitude1{IEXP}{1} - amplitude2{IEXP}{1})./amplitude1{IEXP}{1};
plot(xplot, yplot, 'ko-')
title('Relative amp diff')
xlim(setaxlim(xplot, 1))

% Phase difference from expected
nexttile(6)
hold on
box on
yplot = (phase1{IEXP}{1} - 180) - (phase2{IEXP}{1} - 0);
plot(xplot, yplot, 'ko-')
title('Phase difference from expected')
xlim(setaxlim(xplot, 1))

%% ------------------------- STANDING ESE --------------------------------
% Baseline correction
OptBl.polyOrder = 2;
OptBl.width = 0.1;
for iexp = 1:nExp
    for ii = 1:nMeas
        [ypbl{iexp}{ii}, blp{iexp}{ii}] = ...
            correctbaseline(xx{iexp}, yp{iexp}{ii}, OptBl);
        [ymbl{iexp}{ii}, blm{iexp}{ii}] = ...
            correctbaseline(xx{iexp}, ym{iexp}{ii}, OptBl);
    end
end

% Integration - Get amplitude and phase
for iexp = 1:nExp
    [valMax, jMax] = max(real(ypbl{iexp}{1}));
    [~, jWidth] = min(abs(real(ypbl{iexp}{1}(jMax:end)) - valMax/10));
    for ii = 1:nMeas
        yp1{iexp}(ii) = trapz(ypbl{iexp}{ii});
        ym1{iexp}(ii) = trapz(ymbl{iexp}{ii});
    end
    ypamp{iexp} = abs(yp1{iexp});
    ypph{iexp} = angle(yp1{iexp})*180/pi;
    ymamp{iexp} = abs(ym1{iexp});
    ymph{iexp} = angle(ym1{iexp})*180/pi;
end

IEXP = 1;
plotblcorrection(6, xx{IEXP}, yp{IEXP}, ym{IEXP}, blp{IEXP}, blm{IEXP});
plotstandingese(12, pulseAmps, ypamp{IEXP}, ypph{IEXP}, ymamp{IEXP}, ymph{IEXP});

%% ------------------------- SAVE TO TXT ---------------------------------
% ------------------------------------------------------------------------
IEXP = 2;
chname = ["mpfu_BrX", "mpfu_BrMinX"];
for ip = 1:2
    savepath = chname(ip) + "_ampPhase_2025-03-28.txt";
    savemat = [mpfuAmps; amplitude{IEXP}{ip}'; phase{IEXP}{ip}'];
    % fid = fopen(savepath, 'w');
    % fprintf(fid, "%f %f %f\n", savemat);
    % fclose(fid);
end

%% ------------------------- FUNCTIONS -----------------------------------
function plotblcorrection(figNum, x, yp, ym, blp, blm)
    nMeas = numel(yp);
    figure(figNum)
    f.Name = "Stanging-ESE: echoes";
    clf
    hold on
    box on
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    cmap = viridis(nMeas);
    funplots = ["real", "imag"];
    for iplot = 1:2
        nexttile(iplot)
        hold on
        box on
        for ii = 1:nMeas
            yplot = feval(funplots(iplot), yp{ii});
            plot(x, yplot, 'Color', cmap(ii, :))
            yplot = feval(funplots(iplot), ym{ii});
            plot(x, yplot, 'Color', cmap(ii, :))
        end
        title(funplots(iplot))
        % xline([x(jMax), x(jMax - jWidth), x(jMax + jWidth)])
        xlim([100, 350])
    end
    for ii = nMeas-6:nMeas
        nexttile(ii - nMeas + 9)
        hold on
        box on
        yplotr = real(yp{ii});
        yploti = imag(yp{ii});
        plot(x, yplotr)
        plot(x, yploti)
        yplotr = real(ym{ii});
        yploti = imag(ym{ii});
        plot(x, yplotr)
        plot(x, yploti)
        yplotblp = imag(blp{ii});
        yplotblm = imag(blm{ii});
        plot(x, yplotblp, 'r')
        plot(x, yplotblm, 'r')
    end
end

function plotstandingese(figNum, x, ypamp, ypph, ymamp, ymph)
    f = figure(figNum);
    f.Name = "Stanging-ESE: Amplitude and phase";
    tL = tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact");
    % Amplitude
    nexttile(1)
    hold on
    box on
    plot(x, ypamp, 'o-')
    plot(x, ymamp, 'o-')
    
    % Relative amplitude diff
    nexttile(3)
    yplot = (ypamp - ymamp)./ypamp;
    plot(x, yplot, 'ko-')
    
    % Phase
    nexttile(2)
    hold on
    box on
    plot(x, ypph, 'o-')
    plot(x, ymph, 'o-')
    
    % Phase diff from expected
    nexttile(4)
    yplot = (ypph - 0) - (ymph + 180);
    plot(x, yplot, 'ko-')
end