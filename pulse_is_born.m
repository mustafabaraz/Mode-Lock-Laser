pt=10000; % No. of time (T) steps.
z_steps=3000; % No. of roundtrips.
window=200; % No. of tau points.
d_T=1/1000; % Time (T) step.
factor=1.763; % Associated with the secant funtion.
T=1000/factor; % To set T in fs.
tau= (-pt/2:pt/2-1); % Tau array in seconds.
f=[(0:pt/2-1) (-pt/2:-1)]; % Frequiency array in 1/s.
w=2.*pi.*f*(window/pt); % Angular frequiency array in rad/s.
l_central = 1550; % Central wavelength (thanks to Prof. Omer Ilday's auxilary code)
fr = (300000/l_central + fftshift(f)*(window/T)*1000.0/pt); % Frequiency in 1/fs.
lambda = (1000*300.0./fr); % Wavelength in nm.
% Modified parameters (different from given in HW4).
beta=1;
D_g=1/500^2;
D_d=0.1;
g=1;
l=0.4;
K=10/3;
delta=2.6;
N=d_T*(K+1i*delta); % Nonlinear operator.
a=0; % Uniform solution= 0.
noise = rand(1,pt); % Fluctuations.
noise_factor = 1E-30; % Fraction of total energy
a = a + noise/sum(noise)*noise_factor; % Fluctuated solution (but not stable).
input_pulse=a;
% Two loops as indicated in HW4.
for n=1:z_steps
 s=n; % To make the a (i.e. input_pulse) array a matrix for mesh plotting.
 for t_r=0:3
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
 input_pulse(s,:)=a; % Fill the rows.
 end
 input_pulse(s+1, :)=input_pulse(s,:); % Update the row values.
end
ttt= ((-pt/2:pt/2-1).*T)./window; % ttt stands for the time array in femtoseconds ( remember that tau was in seconds)
% 3-D waterfall plot.
mesh(ttt, 1:z_steps+1, abs(input_pulse).^2)
xlabel('Time (fs)')
ylabel('No. of Roundtrips')
zlabel('Intensity (a.u.)')
