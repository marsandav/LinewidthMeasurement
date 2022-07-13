%% Read ITF

[file,folder] = uigetfile('*.txt','Select ITF file');
fn=[folder,file];
ITF=load(fn); 

power = 0.2;                    % Optical power used to measure the ITF [mW]
Gdc2ac = 9;                     % Gain conversion DC to AC of photodiode
dVdP=max(ITF(:,2))/power;       % Conversion factor from intensity to voltage


% ITF smoothing (prevents over-estimation of dIr)
sf = 11;                        % Smoothing factor
ITF(:,2)=gaussianSmoothing1D(ITF(:,2),sf);

% First derivative of ITF
dIr=diff((ITF(:,2)./max(ITF(:,2))))./diff(ITF(:,1));
[maxdIr,~]=max(abs(dIr(:)));

%% Frequency noise calculation

folder=uigetdir('C:\');
files=ls([folder,'\*.txt']);

P0=14.3;                        % Optical power of the laser output
termination = 2;                % Factor compensation 50 ohm termination
c=299792458;                    % Light speed [m/s]
lamb_centre=1550e-9;            % Central wavelength [m]
f_centre=c/lamb_centre;         % Central frequency [Hz]
    
k=1;
for i=1:numel(files(:,1))
    fn=[folder,'\',files(i,:)];
    D=load(fn);
    
    time=(D(:,1)+abs(D(1,1)));  % Normalise time scale and convert to [s] 
    waveform=D(:,2);            % Noise waveform in [V]    

    dl = 1e-9*termination*waveform./(maxdIr*dVdP*Gdc2ac*P0);    % Wavelength fluctuation in [m]
    df = dl*c/(lamb_centre^2);                                  % Frequency fluctuation in [Hz]

    % Fourier transform
    dt = time(2);
    fs = 1/dt;
    N = numel(time);
    xdft = fft(df);
    xdft = xdft(1:N/2+1);
    psdx = (1/(fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fs/length(D*1e3):fs/2;
    PSD(k,:)=sqrt(psdx);

    k=k+1;
end

%% Plot averaged sqrt PSD of frequency noise
figure('color','w','DefaultAxesFontSize',14);
loglog(freq,mean(PSD,1))
grid on
xlabel('Frequency (Hz)')
ylabel('Frequency noise (Hz/sqrt(Hz))')     
