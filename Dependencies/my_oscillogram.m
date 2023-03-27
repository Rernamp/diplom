function no = my_oscillogram(IQ_mas, SAMPLE_RATE_Hz)

len = length(IQ_mas);
t_mas = (0:1:len-1)/SAMPLE_RATE_Hz*1000;
figure();
hold all;
plot(t_mas, real(IQ_mas));
plot(t_mas, imag(IQ_mas));
plot(t_mas, abs(IQ_mas));
hold off;
xlabel('t, мсек');
title('Осциллограмма');
grid on;

no = 1;

end

