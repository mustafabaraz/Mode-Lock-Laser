MATLAB Code:
pt=10000; % No. of time (T) steps.
z_steps=1000; % No. of roundtrips.
window=150; % No. of tau points.
d_T=1/100000; % Time (T) step.
factor=1.763; % Associated with the secant funtion.
T=1000/factor; % To set T in fs.
tau= (-pt/2:pt/2-1); % Tau array in seconds.
f=[(0:pt/2-1) (-pt/2:-1)]; % Frequiency array in 1/s.
w=2.*pi.*f*(window/pt); % Angular frequiency array in rad/s.
l_central = 1550; % Central wavelength (thanks to Prof. Omer Ilday's auxilary code)
fr = (300000/l_central + fftshift(f)*(window/T)*1000.0/pt); % Frequiency in 1/fs.
lambda = (1000*300.0./fr); % Wavelength in nm.
% Parameters given in HW 4.
beta=0.95;
D_g=1/500^2;
D_d=-1;
g=1;
l=0.4;
K=2/3;
delta=2;
N=d_T*(K+1i*delta); % Nonlinear operator.
a=sech(tau./window).^(1+1i*beta); % Initial pulse.
% Two loops as indicated in HW4.
for n=1:z_steps
 for t_r=0:10
 % First take a half step with dispersion only.
 a_f=fft(a);
 a_f=a_f.*exp((d_T/2)*((g-l)+(D_g-1i*D_d)*-(w).^2));
 a=ifft(a_f);

 % Then take a full step with nonlinear operator only.
 a=a.*exp(N.*(a.*conj(a)));

 % Finally take a half step with dispersion only again.
 a_f=fft(a);
 a_f=a_f.*exp((d_T/2)*((g-l)+(D_g-1i*D_d)*-(w).^2));
 a=ifft(a_f);
 end

end
% Plot the initial and final pulse shapes.
subplot(3,1,1)
hold on
ttt= ((-pt/2:pt/2-1).*T)./window; % ttt stands for the time array in femtoseconds ( remember that tau was in seconds)
input_pulse=sech(tau./window).^(1+1i*beta); % define the initial pulse as in Eq. (31).
plot(ttt, abs(a).^2, '--k')
axis([-1000 1000 0 inf]);
plot(ttt, abs(input_pulse).^2, '-r')
title('Shapes of the pulses inside the cavity')
xlabel('Time (fs)')
ylabel('Intensity')
legend('Shape of the initial pulse.','Shape of the final pulse.')
hold off
% Plot the power spectra of initial and final pulses.
subplot(3,1,2)
hold on
plot(lambda,abs(fftshift(fft(a))).^2./max(abs(fftshift(fft(a))).^2), '--k') %
Notice the normalization.
axis([1520 1580 0 inf]);
plot(lambda,abs(fftshift(fft(input_pulse))).^2./max(abs(fftshift(fft(input_pulse))).^2), '-r') % Notice the normalization.
title('Power spectra of the pulses')
xlabel('Wavelength (nm)')
ylabel('Intensity')
legend('Spectrum of the initial pulse','Spectrum of the final pulse')
hold off
