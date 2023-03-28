function taps = fir_rcos(Interp_factor, N_slow_samples_one_side, rolloff)

   taps = rcosdesign(rolloff, 2*N_slow_samples_one_side, Interp_factor, 'normal');

%   taps = rcosdesign(rolloff, 2*N_slow_samples_one_side, Interp_factor, 'sqrt');

   % Нормировка, чтобы усиление на частоте 0 Гц было 0 дБ:
   taps = taps / sum(taps);

end % function
