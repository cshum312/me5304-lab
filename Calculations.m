data = readtable('datamatlab.csv','NumHeaderLines',1);

global A B x4 x3 x2 x1 x0 n p u_inf
A = 2.1941;
B = 0.7496;

x4 = 159.41;
x3 = -1178.4;
x2 = 3276.3;
x1 = -4037;
x0 = 1853;

p = [x4 x3 x2 x1 x0];

n = 132000;

cali = table2array(data(:,"Var2"));
sample6khz = table2array(data(:,"Var3"));
sample300hz = table2array(data(:,"Var4"));
sample400hz = table2array(data(:,"Var5"));

u_inf = 12.6;

sample300hz = rmmissing(sample300hz);
sample400hz = rmmissing(sample400hz);

figure(1)
[mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(cali, 'poly', 6000);
sgtitle('Calibration - Polynomial')

figure(2)
[mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(cali, 'king', 6000);
sgtitle('Calibration - Kings Law')

figure(3)
[mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(sample6khz, 'poly', 6000);
sgtitle('6000Hz - Polynomial')

figure(4)
[mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(sample6khz, 'king', 6000);
sgtitle('6000Hz - Kings Law')

figure(5)
[mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(sample300hz, 'poly', 300);
sgtitle('300Hz - Polynomial')

figure(6)
[mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(sample300hz, 'king', 300);
sgtitle('300Hz - Kings Law')

figure(7)
[mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(sample400hz, 'poly', 400);
sgtitle('400Hz - Polynomial')

figure(8)
[mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(sample400hz, 'king', 400);
sgtitle('400Hz - Kings Law')


function [mean_val, var_val, rms_val, tu_val, skew_val, kurt_val] = process(data, type, sampling)
    global p A B u_inf

    if type == 'poly'
        U = polyval(p,data);
    elseif type == 'king'
        U = ((1/B)*(data.^2 - A)).^(1/0.45);
    else
        return
    end

    xStart = 0;
    dx = 1/sampling;

    if sampling == 6000
        N = 132000;
        t = xStart + (0:N-1)*dx;
    else
        N = 10000;
        t = xStart + (0:N-1)*dx;
    end

    mean_val = mean(U);
    var_val = var(U);
    rms_val = std(U);
    tu_val = rms_val/mean_val;
    skew_val = skewness(U);
    kurt_val = kurtosis(U);
    
    subplot(2,2,1)
    histogram(U, 'Normalization', 'pdf');
    hold on;
    ax = gca;
    x = linspace(ax.XLim(1), ax.XLim(2), 1000);
    plot(x, pdf(makedist('Normal', 'mu', mean_val, 'sigma', rms_val), x), 'LineWidth', 2)
    xlabel('Velocity')
    ylabel('PDF')
    title('PDF')

    [acf,lags] = autocorr(U, NumLags=30);
    [acfP,lagsP] = autocorr(U, NumLags=1000);
    temp = t(1:length(lags));
    tempP = t(1:length(lagsP));
    
    subplot(2,2,2)
    plot(temp, acf);
    xlabel('Delay')
    ylabel('Autocorrelation')
    title('Autocorrelation')

    Fs = sampling;
    L = (t(end) + dx)*Fs;
    Y = fft(U);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(5:end-1) = 2*P1(5:end-1);
    freq = Fs*(0:(L/2))/L;

    subplot(2,2,3)
    plot(freq(5:end),P1(5:end)) 
    xlabel('frequency (Hz)')
    title('Vortex Shedding')

    pks = max(P1(5:end));
    f = freq(P1 == pks);
    period = 1/f;

    FsP = sampling;
    LP = (tempP(end) + dx)*FsP;
    YP = fft(acfP);
    P2P = abs(YP/LP);
    P1P = P2P(1:LP/2+1);
    P1P(2:end-1) = 2*P1P(2:end-1);
    freqP = FsP*(0:(LP/2))/LP;

    subplot(2,2,4)
    plot(freqP(5:end),P1P(5:end)) 
    xlabel('frequency (Hz)')
    title('Power Spectrum')

    pksP = max(P1P(5:end));
    fP = freqP(P1P == pksP);
    St = 0.010*fP/u_inf;

end


