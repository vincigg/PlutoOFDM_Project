%% Clearing variables
clear;close all;clc

%% OFDM MODULATOR
% Root of raised cosine filter parameters
alpha = 0.25;    % Rolloff factor
L     = 30;      % Truncated impulse length (-L/2:L/2)
Nc    = 10;      % Number of samples per symbol

% TX filter design
istrt = floor(L/2); 
n     = -istrt:1/Nc:istrt; %(-L/2:1/Nc:L/2)
pt    = truncRRC(n, alpha, 0);

% RX filter design: p(t) = p(-t)
pmt = pt;

% Amplitudes of the paths
a = [0.8; 0.4; 0.4];

% Delays of the paths
tau = [0; -1; 2];

% Frequency-dispersive channel
ptau = sum(a.*truncRRC(n, alpha, tau), 1)/(sum(a));

% Computation of the number of useful subcarriers
N   = 256;                   % IFFT order (number of subcarriers)
Na  = floor(N/2*(1-alpha));  
Nu  = 2*Na + 1;             % Number of useful subcarriers
Nsc = N - Nu;               % Number of suppressed carriers
Np  = 30;                   % Cyclic prefix length


% M-PSK consellation
M      = 4;                         % Number of PSK levels
mb     = log2(M);                   % Number of bits
psklev = exp(1j*2*pi*(0:M-1)/M);    % M-PSK levels
nSym   = 1e3;                      % Number of symbols

% Generation of data input and symbol mapping
symtx_data       = randi([0,M-1], nSym, 1); % Generate random M-ary symbols

% Generation of sequence of data used for equalization

training_equalizer = randi([0,M-1], Nu, 1);

symtx = [training_equalizer; symtx_data];
[symtxg, ~] = togray(symtx, mb);      % Gray mapping
ci          = psklev(symtxg+1);       % Symbol selection

% Grey mapping only training sequence to obtain known ci symbols at
% receiver side to use for equalization
[tr_eq_gray, ~] = togray(training_equalizer,mb);
known_ci = psklev(tr_eq_gray+1);

% Set the length of the PN sequence used for timing estimation
pn_length = N/2;

% Generate a random binary sequence
pn_sequence = -1 + 2.*randi([0, 1], 1, pn_length);

% Repeat the PN sequence to create the training symbol
training_symbol = repmat(pn_sequence, 1, 2);

% Add cyclic prefix to the training symbol
training_symbol_with_prefix = [training_symbol(end - Np + 1:end), training_symbol];

% Serial to parallel block
r = rem(nSym, Nu);
if r~=0
    txt = sprintf('discarding last %d symbols', r);
    disp(txt);
    ci = ci(1:end-r);
end

ci_par = reshape(ci, Nu, []);

Nbl = size(ci_par,2); % Overall number of blocks

% Virtual carrier insertion
ci_N = [ci_par(1:Na+1,:); zeros(Nsc, Nbl); ci_par(Na+2:end,:)];

% IDFT
a_N = N*ifft(ci_N);

% Cyclic prefix insertion
a_NT = [a_N(end-Np+1:end,:); a_N];

% Parallel to series conversion
a_NT = a_NT(:);

% Metric training symbol insertion in time domain
a_NT = [training_symbol_with_prefix.'; a_NT];

% TX filtering: upsample and filter
a_up    = upsample(a_NT, Nc)/Nc;
txSig    = filter(ptau, 1, a_up);


%% Additional trasmission noise

  %SNR     = 10;          % SNR in dB
  %SNRlin  = 10^(SNR/10); % SNR linear scale
  %swn     = sqrt(0.5*Nc/(SNRlin*mb)); % Noise variance
  
  %txSig = txSig + swn*(randn(size(txSig)) + 1j*randn(size(txSig)));

%-------------------------------------------------------------------------
%% Radio paramaters

% Overflow only in SA
SampleRate = 1e6;
SamplesPerRXFrame = 2^14; % Change to 2^16 to remove overflow
FramesToCollect = 5;
FrequencyOffset = 0;

%% Set up radio and capture some data
rx = sdrrx('Pluto','SamplesPerFrame',SamplesPerRXFrame,...
    'BasebandSampleRate',SampleRate,'OutputDataType','double');
tx = sdrtx('Pluto','Gain',-30,...
    'BasebandSampleRate',SampleRate);
tx.CenterFrequency = tx.CenterFrequency + FrequencyOffset;
tx.transmitRepeat(txSig);
% Get data from radio
saved = zeros(FramesToCollect*SamplesPerRXFrame,1);
ofe = 0;
for g=1:SamplesPerRXFrame:FramesToCollect*SamplesPerRXFrame
    [data1,len,of] = rx();
    saved(g:g+SamplesPerRXFrame-1) = data1;
    if of % Count overflows
        ofe = ofe + 1;
    end
end
fprintf('Overflow events: %d of %d\n',ofe,FramesToCollect)

rxSig = saved;

%-------------------------------------------------------------------------

%% Receiver Side
rxdem   = filter(pmt, 1, rxSig);
rx_ds   = downsample(rxdem, Nc);   % downsample demodulated signal
rx_sym = rx_ds(L+1:end);               % remove filter transient

%% Timing Sync Method A

Ltrain = N/2;

finestra = 4000; % vector containing time indexes of the metric. it's empiric, this value could be calculated though

Pgrande = zeros(1, finestra);%empty vectors
Rf = zeros(1, finestra);
Mf = zeros(1, finestra);
M1 = zeros(1, finestra); %metric vector

 for d = 1: finestra
    r_asterisco = conj(rx_sym(d:d+Ltrain-1));
    r_base = rx_sym(d+Ltrain:d+2*Ltrain-1);
    Pgrande(d)= sum(r_asterisco .* r_base);

    % r(d+m) Definition when the sum reaches N-1
    r_N = rx_sym(d:d+N-1);
    Rf(d) = 0.5*sum(r_N.*conj(r_N));


    Mf(d) = (abs(Pgrande(d))^2) / (Rf(d)^2);

 end
 risultato = zeros(1, finestra);

 for d= 1:finestra
    for k = 0 : Np
        if d-k >= 1
            risultato(d) = risultato(d) + Mf(d-k);
        end
    end
 end

M1 = (1/(Np+1)) .* risultato; 


figure, plot(M1), title("M1, Frame syncronization metric")


%frame recognition and extraction

[picco,idx]=max(M1);
[piccok,idxk]=maxk(M1,50);

for indicemassimo = 1:50
    if abs(idxk(indicemassimo)-idx)>300
        secondomax = idxk(indicemassimo);
        break;
    end
end



if idx >secondomax
    change = idx;
    idx = secondomax;
    secondomax=change;
end

primomax = idx;


% Selecting only the correct part of the received signal, corrisponding to
% the frame we transmitted

rx_sym = rx_sym(primomax+N:secondomax-Np-1);

%% S/P and removing CPs
r = rem(length(rx_sym), N+Np);
a_NTr = rx_sym(1:end-r);
a_NTr = reshape(a_NTr, N+Np, []);

a_Nr  = a_NTr(Np+1:end,:); % Removing all Cyclic prefixes

%% FFT
ci_Nr = 1/N*(fft(a_Nr));

%% Removing non useful subcarriers 
ci_parr = [ci_Nr(1:Na+1,:); ci_Nr(end-Na+1:end,:)];

%% Equalization (assuming known channel)
eq_length = Nu; %length of the known symbols used for channel eq

Hstima = [ci_parr(:,1)].' ./ known_ci;

Nblock_eq = size(ci_parr,2); % Overall number of blocks to equalize

equalized = zeros(Nu, Nblock_eq-1);

sym_da_eq_par = ci_parr(:,2:end);

equalized = sym_da_eq_par ./ Hstima.';



ci_eqserie = equalized(:); %P2S conversion
cis = ci_parr(:); %non equalized p2s for reference

figure, scatter(real(ci_eqserie), imag(ci_eqserie)), title("Equalized Symbols, (removed PN sequence)")

figure, scatter(real(cis), imag(cis)), title("NON Equalized Symbols")


%% DETECTION Algorithm
% ML detection (minumum distance)
dist_vec     = abs(psklev(:) - ci_eqserie.').^2;
[~, sym_idx] = min(dist_vec);
det_symg     = sym_idx - 1; % Gray-coded detected symbol
[det_sym, ~] = fromgray(det_symg, mb); % Detected symbol


%% Error valutation
err = symtx_data(1:size(det_sym,1)) == det_sym;
correct_detection = sum(err);
bad_detect = size(det_sym,1) - correct_detection;
error_rate = bad_detect/size(det_sym,1);
fprintf('Error rate: %.2f percent \n', error_rate*100);

%% PAPR
papr_Test = calculatePAPR(rxSig);

% Display PAPR
fprintf('PAPR: %.2f dB\n', papr_Test);

%% IQ Compontents
iq(1,:) = real(txSig);
iq(2,:) = imag(txSig);
% Create a plot of the I-Q components
figure,
plot(iq(1,:),iq(2,:), "-"),title("IQ Diagram")

%% Power Spectrum
% Calculate power spectrum using pwelch, txsig before channel insertion but after
% oversampling

txSigNoF = filter(pt, 1, a_up); %Without frequency selectivity (before encountering the channel)
[pxx, f] = pwelch(txSigNoF,[],[],[],SampleRate);

% Plot the power spectrum
pxx_shifted = circshift(pxx, size(f,1)/2);
fplot = f- ones(size(f,1),1).*max(f)/2;
fplot = fplot/SampleRate;
figure, plot(fplot, 10*log10(pxx_shifted)), xlabel('Normalized Frequency f/Fs'), ylabel("PSD (dB)");
xticks(-0.5:0.1:0.5);  % Set x-axis tick values from -0.5 to 0.5 with a step of 0.1

title('Power Spectrum'), grid on;

%% Functions definition
function papr = calculatePAPR(signal)
    % Calcola il PAPR di un segnale
    % Calcola la potenza di picco
    peakPower = max(abs(signal))^2;
    % Calcola la potenza media
    averagePower = mean(abs(signal))^2;
    % Calcola il PAPR
    papr = 10 * log10(peakPower / averagePower);
end

function pRRC=truncRRC(x,a,tau)
%truncRRC: Truncated raised root cosine impulse
%Use:	pRRC=truncRRC(x,a,tau)
%		x: normalised time
%		a: rolloff, within [0,1]
%       tau: input delay
%		truncRRC(x,a)= (1-a) sinc( (1-a) x) +
%		a( sinc(ax + 1/4) cos(pi (x+1/4)) + sinc(ax - 1/4) cos(pi (x-1/4))

x = x - tau; % delayed time
pRRC = (1-a)*sinc(x*(1-a));
pRRC = pRRC + a* sinc(a*x + 0.25).*cos(pi * (x + 0.25));
pRRC = pRRC + a* sinc(a*x - 0.25).*cos(pi * (x - 0.25));
end

function [gd, gb] = togray(in, len)

if ischar(in)
in = bin2dec(in);
end

gd = bitxor(in, bitshift(in,-1));
gb = dec2bin(gd, len);

end

function [de, bi] = fromgray(in, len)

if ischar(in)
in = bin2dec(in);
end

in = in(:);
n_syms = length(in);

b = zeros(n_syms,len);

% get MSB of in
b(:,1) = bitget(in, repmat(len,[n_syms,1]));

for sh = 2 : len
    in_loc = bitget(in, repmat(len-sh+1,[n_syms,1]));
    b(:,sh) = bitxor(b(:,sh-1), in_loc);
end

bi = num2str(b);
de = bin2dec(bi);
end