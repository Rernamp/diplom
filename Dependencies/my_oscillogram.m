function no = my_oscillogram(IQ_mas, SAMPLE_RATE_Hz)

len = length(IQ_mas);
t_mas = (0:1:len-1)/SAMPLE_RATE_Hz*1000;
figure();
hold all;
plot(IQ_mas);
hold off;
grid on;

no = 1;

end

