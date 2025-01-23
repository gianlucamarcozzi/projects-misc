clearvars

squarePulse = @(xx, p) (p(1) < xx & p(2) > xx).*(1 - exp(-(xx - p(1))/500)) + ...
    (~(p(1) < xx | p(2) > xx))*0;

nx = 10001;
xx = linspace(0, nx - 1, nx);
p1 = [2000, 3000];
yy = squarePulse(xx, p1);
p2 = [2000, 3100];
yy2 = squarePulse(xx, p2);
% yy = yy*800;

figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
nexttile
plot(xx, yy)
hold on
plot(xx, yy2)

% FFT
tStep = xx(2) - xx(1);
fSampl = 1/tStep;
nzf = 0;  % Zero filling
if nzf <= nx && nzf ~= 0
    warning("nzf <= nx. Continuing without zero-filling.")
    nzf = 0;
end
if nzf == 0
    ff = fSampl/(nx)*(-(nx)/2:(nx)/2 - 1);
else
    ff = fSampl/(nzf)*(-(nzf)/2:(nzf)/2 - 1);
end


% tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")


if nzf ~= 0  % Zero filling
    yy(nzf) = 0;  % Zero filling
    yy2(nzf) = 0;  % Zero filling
else  % No zero filling (adjust values in memory in case)
    yy = yy(1:nx);
    yy2 = yy2(1:nx);
end

fp0 = fft(yy);
fp2 = fft(yy2);

% fp0(2) = fp0(2)*2/3;
% fp0(end) = fp0(end)*2/3;
% fp0(1) = fp0(1)*3/4;
nexttile(2)
plot(ff, abs(fftshift(fp0)), 'o-')
hold on
plot(ff, abs(fftshift(fp2)), 'o-')

xlim(0.003*[-1, 1]*1)

yb = ifft(fp0);
nexttile(1)
hold on
plot(xx, yb)
