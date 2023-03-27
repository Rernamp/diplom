subsystemType = '32APSK 3/4';

addpath("Dependencies");

fs = 1;

EsNo = 20;

numSymbols = 1e4;

a_A = 2.1587; 
b_A = 1.1517;
a_F = 4.0033;
b_F = 9.1041;

dvb = getParamsDVBS2Demo(subsystemType, EsNo, 50);


sps = 8; % samples per symbol
SAMPLE_RATE_Hz = 48000;
Time_sec = 5;
bps = dvb.BitsPerSymbol;
N_symbols = Time_sec * SAMPLE_RATE_Hz / sps;
N_bits = N_symbols * bps;



pskModulator = comm.PSKModulator(...
'ModulationOrder', dvb.ModulationOrder,...
'BitInput', true, ...
'PhaseOffset', dvb.PhaseOffset, ...
'SymbolMapping', 'Custom', ...
'CustomSymbolMapping', dvb.SymbolMapping);

if dvb.ModulationOrder == 4 || dvb.ModulationOrder == 8
const = pskModulator.constellation;
else
const = dvbsapskmod((0:dvb.ModulationOrder-1)', dvb.ModulationOrder, 's2', ...
  dvb.CodeRate, 'UnitAveragePower', true);
end
maxConstReal = max(real(const));
maxConstImag = max(imag(const));

sigmas = 10^(-EsNo/10);
maxReal = maxConstReal + 3*sqrt(sigmas);
minReal = -maxReal;
maxImag = maxConstImag + 3*sqrt(sigmas);
minImag = -maxImag;

% 1.1) Передаваемое сообщение:
mas_Tx_message = randi(2, N_bits, 1)-1;

% 1.2) Формирвание символов из бит: 
matr_Tx_message = reshape(mas_Tx_message, N_symbols, bps);  % Reshape data into binary n-tuples
mas_Tx_int_symbols = bi2de(matr_Tx_message);                                     % Convert to integers

% 1.3) Формирование символов созвездия:
mas_Tx_clx_symbols = zeros(N_symbols, 1);
for i = 1 : 1 : N_symbols
   mas_Tx_clx_symbols(i) = const(mas_Tx_int_symbols(i)+1);
end % i


xComplex = mas_Tx_clx_symbols;
% xComplex = (rand(numSymbols,1)*(maxReal-minReal)+minReal) + ...
% 1i*(rand(numSymbols,1)*(maxImag-minImag)+minImag);

% sig_in = abs(xComplex) .* exp(1i * 2 * pi * angle(xComplex));
sig_in = xComplex;
scatterplot(sig_in)

rolloff = 0.5;
FIR_h = fir_rcos(sps, 3, rolloff); % формирующий фильтр интерполятора

sig_in  = sps * conv(sig_in, FIR_h, 'same');

my_oscillogram(sig_in, 1);


G = a_A.*abs(xComplex)./(1+b_A.*(abs(xComplex).^2));
F = a_F*(abs(xComplex).^2)./(1+b_F*(abs(xComplex).^2));
sig_out1 = G .* exp(1i * 2 * pi .* F);
sig_out = sig_in .* (sig_out ./ abs(sig_in));
%%sig_out_Sale = G.*exp(1i * 2 * pi .* (F + angle(xComplex)));

my_oscillogram(sig_out, 1);

peakFactor(xComplex)
spec_dB(sig_in, 1, "Input") 
spec_dB(sig_out, 1, "Output Saleh 32APSK")
%%spec_dB(sig_out_Sale, 1, "From Saleh sci-hub")

load('h_FIR_Rx.mat');

mas_Rx_IQ = conv(sig_out, h_FIR_Rx, 'same');
mas_Rx_clx_symbols = mas_Rx_IQ(1 : sps : end);

mas_Rx_int_symbols = zeros(N_symbols, 1);
for i = 1 : 1 : N_symbols
  Rx_clx_symbol = mas_Rx_clx_symbols(i);
  
  [vmin imin] = min(abs(const - Rx_clx_symbol));
  
  mas_Rx_int_symbols(i) = imin-1;
end % for i

matr_Rx_message = de2bi(mas_Rx_int_symbols, bps);
mas_Rx_message = matr_Rx_message(:);

my_oscillogram(mas_Rx_message, 1);



function [peakFactor] = peakFactor( xComplex)
%     peakFactor = db((max(signal.^2)) ./ (mean(signal.^2)));
    peakFactor = db(peak2rms(xComplex))
end

function [ spec ] = spec_dB( sig_out , fs, titleName)

    x_fft = fft(sig_out);
    
    ahc_x = fftshift(mag2db(abs(x_fft)));
    
    f = (0:length(sig_out)-1)*fs/length(sig_out);
    
    figure();
    plot(f, ahc_x);
    set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
    set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman'); 
    grid on;
    title(titleName);
    xlabel('Частота, Hz');
    ylabel("BP, дБ")
end





