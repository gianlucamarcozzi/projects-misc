clearvars

%% IMPORT

generalFolder = "../data/raw/";
expName = "misc001-E-002";
measFolder = generalFolder + expName;

[x00, y0, Param0] = loadfolderelexsys(measFolder, '*1.DTA');
nMeas = numel(Param0);
x0 = x00{1};
[x22, y2, Param2] = loadfolderelexsys(measFolder, '*2.DTA');
x2 = x22{1};

temp = readtable(measFolder + "/" + expName + "-param.txt");
pulseAmp = temp{:, 2};

%% ECHO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normModel = @(xx, p) p(1)*exp(-(xx - p(2)).^2/2/p(3)^2);
fitmodel = @(p) normModel(x0, p);
p0 = [0, 250, 40];

% Smooth
for ii = 1:nMeas
    % y{ii} = datasmooth(y0{ii}, 10, 'savgol');
    y{ii} = y0{ii};
end
fitOpt = optimoptions('lsqnonlin','Display','off');
for ii = 1:nMeas
    % Real
    ydata = real(y{ii});
    p0(1) = max(abs(ydata));
    if -min(ydata) > max(ydata)
        p0(1) = -p0(1);
    end
    [pfit{ii, 1}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pci(ii, 1, :, :) = ...
        nlparci(pfit{ii, 1}, residual, 'jacobian', jacobian);
    [yfit{ii}, pfit{ii, 1}(1), pfit{ii, 1}(4)] = mldividefun(...
        fitmodel, ydata, pfit{ii, 1});

    % Imag
    ydata = imag(y{ii});
    p0(1) = max(abs(ydata));
    if -min(ydata) > max(ydata)
        p0(1) = -p0(1);
    end
    [pfit{ii, 2}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pci(ii, 2, :, :) = ...
        nlparci(pfit{ii, 2}, residual, 'jacobian', jacobian);
    [yfitImag, pfit{ii, 2}(1), pfit{ii, 2}(4)] = mldividefun(...
        fitmodel, ydata, pfit{ii, 2});
    yfit{ii} = yfit{ii} + 1i*yfitImag;
end

for ii = 1:nMeas
    % vre = max(real(yfit{ii}));
    % vim = min(imag(yfit{ii}));
    vre(ii) = trapz(real(yfit{ii}) - pfit{ii, 1}(4));
    vim(ii) = trapz(imag(yfit{ii} - 1i*pfit{ii, 2}(4)));
    
    phase(ii) = atan(vim(ii)/vre(ii));
    magnitude(ii) = sqrt(vre(ii)^2 + vim(ii)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yshift = -2000;
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact") 
for ii = 32:51
    nexttile
    plot(x0, real(y{ii}), 'o-')
    hold on
    plot(x0, real(yfit{ii}))
    plot(x0, imag(y{ii}) + yshift, 'o-')
    plot(x0, imag(yfit{ii}) + yshift)
end
% savefigas(gcf, "../images/misc001-E-002_01_exampleFitEchoes.png")

% idxs = [1, 4:13, 15:51];
figure(2)
clf
plot(pulseAmp, phase*180/pi, 'o-')
hold on
% plot(pulseAmp(dxs), phase(idxs), 'o-')
yyaxis right
plot(pulseAmp, magnitude, 'o-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% SAVE TO FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen('output.txt', 'w');
for ii = 1:numel(phase)
    fprintf(fileID, '%.2f, ', phase(ii));
end
fclose(fileID);

%% PHASE CORRECTION

phaseshift = @(y, p) y*exp(1i*p);
for ii = 1:nMeas
    ys{ii} = phaseshift(y{ii}, -phase(ii));
end

fitOpt = optimoptions('lsqnonlin','Display','off');
for ii = 1:nMeas
    % Real
    ydata = real(y{ii});
    p0(1) = max(abs(ydata));
    if -min(ydata) > max(ydata)
        p0(1) = -p0(1);
    end
    [pfits{ii, 1}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pci(ii, 1, :, :) = ...
        nlparci(pfits{ii, 1}, residual, 'jacobian', jacobian);
    [ysfit{ii}, pfits{ii, 1}(1), pfits{ii, 1}(4)] = mldividefun(...
        fitmodel, ydata, pfits{ii, 1});
end

for ii = 1:nMeas
    % vre2(ii) = max(real(ysfit{ii}));
    vre2(ii) = trapz(real(ysfit{ii}) - pfits{ii, 1}(4));
    
    magnitudes(ii) = vre2(ii);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yshift = -2000;
figure(3)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact") 
for ii = 32:51
    nexttile
    plot(x0, real(ys{ii}), 'o-')
    hold on
    plot(x0, real(ysfit{ii}))  % Fits are not perfect at large amplitudes
    plot(x0, imag(ys{ii}) + yshift, 'o-')
    % plot(x0, imag(yfit{ii}) + yshift)
end
savefigas(gcf, "../images/misc001-E-002_04_echoesPhaseCorrected.png")

figure(4)
clf
plot(pulseAmp, magnitude, 'o-')
yyaxis right
plot(pulseAmp, magnitudes, 'o-') % Fits are not perfect at large amplitudes

%% RABI FIT cosine

% rabii = rabii(:, 1:nTau);  % Recover rabii before zero filling
expcos = @(xx, A, tau, w) A*exp(-xx/tau).*cos(2*pi*w*xx);
fitmodel = @(p) expcos(x2, 1, p(1), p(2));

fitOpt = optimoptions('lsqnonlin','Display','off');

w0 = (linspace(sqrt(3), sqrt(48), nMeas)).^2*1e-3;
for ii = 1:nMeas
    % Real
    p0(2) = w0(ii);
    if ii < 17
        p0(1) = 350;
    else
        p0(1) = 400;
    end

    ydata = real(y2{ii});
    
    % Don't fit amplitude in the lsqnonlin
    [pfitrr{ii}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pcirr(ii, :, :) = nlparci(pfitrr{ii}, residual, 'jacobian', jacobian);

    % Get proper yfit, amplitude and constant
    [yfitrr{ii}, pfitrr{ii}(3), pfitrr{ii}(4)] = mldividefun(...
        fitmodel, ydata, pfitrr{ii});
    
    % Imag
    p0(2) = w0(ii);
    if ii < 19
        p0(1) = 250;
    else
        p0(1) = 400;
    end

    ydata = imag(y2{ii});
    
    % Don't fit amplitude in the lsqnonlin
    [pfitri{ii}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pciri(ii, :, :) = nlparci(pfitri{ii}, residual, 'jacobian', jacobian);

    % Get proper yfit, amplitude and constant
    [yfitri{ii}, pfitri{ii}(3), pfitri{ii}(4)] = mldividefun(...
        fitmodel, ydata, pfitri{ii});
end
for ii = 1:nMeas
    bestPi(ii) = 1/2/pfitrr{ii}(2);
    bestPiIm(ii) = 1/2/pfitri{ii}(2);
    experPi(ii) = eval(Param0{ii}.PlsSPELGlbTxt(192:193));

    phaseRabi(ii) = atan(pfitri{ii}(3)/pfitrr{ii}(3));
    magnitudeRabi(ii) = hypot(pfitri{ii}(3), pfitrr{ii}(3));
    
end

for ii = 1:nMeas
    freqre(ii) = pfitrr{ii}(2)*1e3;
    diffreim(ii) = (pfitrr{ii}(2) - pfitri{ii}(2))/pfitrr{ii}(2);

    fprintf("Freq Re: %f, Re - Im: %f %%\n", freqre(ii), diffreim(ii)*100)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yshift = -3e4;
figure(6)
clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
for ii = 32:nMeas
    nexttile
    plot(x2, real(y2{ii}), 'o-')
    hold on
    plot(x2, yfitrr{ii})
    plot(x2, imag(y2{ii}) + yshift, 'o-')
    plot(x2, yfitri{ii} + yshift)
    title(sprintf('%d: %.4f MHz, %.2f ns', ...
          ii, pfitri{ii}(2)*1e3, 1/2/pfitri{ii}(2)))
    plotText = sprintf("%.2f", pulseAmp(ii));
    text(gca, 0.8, 0.8, plotText{end}, 'Units', 'normalized')
    % xline(32)
end
savefigas(gcf, "../images/misc001-E-002_02_exampleFitRabi.png")

figure(7)
clf
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact')
ax1a = nexttile;
plot(pulseAmp, phase*180/pi, 'o-')
hold on
plot(pulseAmp, phaseRabi*180/pi, 'o-')
legend("Phase echo", "Phase Rabi")
% plot(pulseAmp(dxs), phase(idxs), 'o-')
ylabel(ax1a, "Phase / deg")
ax2a = nexttile;
plot(pulseAmp, magnitude, 'o-')
yyaxis right
plot(pulseAmp, magnitudeRabi, 'o-')
xline(79.6)
legend("Magnitude echo", "Magnitude Rabi")
ax3a = nexttile;
plot(pulseAmp, experPi, '-o')
hold on
plot(pulseAmp, bestPi, '-o')
plot(pulseAmp, bestPiIm, '-o')
legend("Experimental pi", "Best pi real(rabi)", "Best pi im(rabi)")
xline(79.6, "HandleVisibility", "off")
xlabel(ax3a, "MPFU yAmp / %")
ylabel(ax3a, "Time length / ns")

% Secondary x-axis
set(ax1a,'box','off')
set(ax2a,'box','off')
set(ax3a,'box','off')
ax1b = axes('Position', ax1a.Position, 'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', 'Color', 'none');
xlim(ax1b, setaxlim(freqre, 1))
xlabel(ax1b, "Rabi freq / MHz")
yticks(ax1b , [])
ax2b = axes('Position', ax2a.Position, 'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', 'Color', 'none');
xlim(ax2b, setaxlim(freqre, 1))
yticks(ax2b , [])
ax2b = axes('Position', ax3a.Position, 'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', 'Color', 'none');
xlim(ax2b, setaxlim(freqre, 1))
yticks(ax2b , [])

savefigas(gcf, "../images/misc001-E-002_03_phase-magnitude-piLength.png")

% xline([-10, 10])

% line(ff2, abs(fftshift(fp0{ii})), 'Parent', ax2(ii), 'Color', 'none');
% 
% ax2(ii).XColor = 'r';
% ax2(ii).XLabel.String = 'Secondary X-axis';
% ax2(ii).XLabel.Color = 'r';

