%%%%%%%%%%%%%%%%%%%%%%%%%%%    1      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read the attached audio file , plot the spectrum of this signal. 
[signal,Fs] = audioread('eric.wav'); % read audio file
signal_samples = length(signal); %number of samples
t= linspace(0,length(signal)/Fs,length(signal)); %generates n points
f= (-Fs/2:Fs/length(signal):Fs/2-Fs/length(signal));% frequency range to zero-center the shifted signal
figure(1)
plot(t,signal)
title ('original signal in time domain')

f_signal = fftshift(fft(signal));
figure (2)
plot(f,abs(f_signal))
title('original signal in freq. domain')


%%%%%%%%%%%%%%%%%%%%%%%%      2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use an ideal filter to remove all frequencies greater than 4KHZ
Fn=Fs/2;
cutoff_frequency = 4000 /Fn; %normalized cutoff frequency
[zeroes, poles] = butter(25, cutoff_frequency, 'low');
t = [0:length(signal_before_modulation)-1]*(1/Fs); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Obtain the filtered signal in both frequency and time domain
filtered_signal_t = filter(zeroes, poles, signal); %Using zeroes and poles of filter with the signal
figure(3);
plot(t, filtered_signal_t);
title('filtered signal in time doain');
xlabel('Time');
ylabel('Amplitude');   


% here filtered signal in freq domain before fftshift
% filtered_signal_f2 = fft(signal_before_modulation); %Change from time domain to frequency domain
% figure(4);
% plot(abs(filtered_signal_f2));
% title('filtered signal in freq domain before fftshift2');
% xlabel('Frequency');
% ylabel('Amplitude');


%f-domain representation after fftshift:
filtered_signal_f =  fftshift(fft(filtered_signal_t), 2);
f = (-Fs/2:Fs/length(signal):Fs/2-Fs/length(signal));   %by default, it plots frequency in bins, that's to change it to Hz 
n = length(signal);
fshift = (-n/2:n/2-1)*(Fs/n); % zero-centered frequency range
figure(5);
plot(fshift, abs(filtered_signal_f));
title('filtered signal in freq domain after fftshift');
xlabel('frequency (Hz)');                                                
ylabel('amplitude');


%%%%%%%%%%%%%%%%%%                4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sound the filtered audio signal

%sound( real(double(filtered_signal_t)), Fs);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modulate the carrier with the filtered signal you obtained
%DSBSC
Fc = 100000;
[N,D] = rat((5*Fc)/Fs);
signal_resampled = resample(filtered_signal_t, N, D);
%resampling after increasing the sampling frequency
time = (0 : 1/(Fs*10) : signal_samples/Fs);
carrierSignal = cos(2*pi*Fc*time);
signal_to_modulate = signal_resampled(1:length(carrierSignal));
DSBSC =  signal_to_modulate.* carrierSignal';

%DSBTC
m=.5;
Mc = max(filtered_signal_t);
normalized_filtered_signal_t=(filtered_signal_t)/Mc;
Ac=2*Mc;
carrierSignal2 = Ac*cos(2*pi*Fc*time);
signal_resampled2 = resample(normalized_filtered_signal_t, N, D);
signal_to_modulate2 = signal_resampled2(1:length(carrierSignal));
DSBTC =  (1+m*signal_to_modulate2).* carrierSignal2';


%DSB-SC ploting
figure(6);
plot(time, DSBSC);                             
title('DSB-SC signal in time domain');
xlabel('time');
ylabel('amplitude');
figure(7);
plot(abs(fft(DSBSC)));                             
title('DSB-SC signal in freq domain');
xlabel('freq');
ylabel('amplitude');

%DSB-TC ploting
figure(8);
plot(time, DSBTC);                             
title('DSB-TC signal in time domain');
xlabel('time');
ylabel('amplitude');
figure(9);
plot(abs(fft(DSBTC)));                            
title('DSB-TC signal in freq domain');
xlabel('freq');
ylabel('amplitude');

% % another solution for step 5
% %DSB_SC
% Fc = 100000;
% [N,D] = rat((5*Fc)/Fs);
% signal_resampled = resample(filtered_signal_t, N, D);
% %resampling after increasing the sampling frequency
% time = (0 : 1/(Fs*10) : signal_samples/Fs);
% carrierSignal = cos(2*pi*Fc*time);
% signal_to_modulate = signal_resampled(1:length(carrierSignal));
% DSBSC =  signal_to_modulate.* carrierSignal';
% figure(6);
% plot(time, DSBSC);                            
% title('DSB-SC signal (resampled)');
% xlabel('time');
% ylabel('amplitude');

% %DSB_TC
% m=.5;
% Mc = max(filtered_signal_t);
% normalized_filtered_signal_t=(filtered_signal_t)/Mc;
% Ac=2*Mc;
% carrierSignal2 = Ac*cos(2*pi*Fc*time);
% signal_resampled2 = resample(normalized_filtered_signal_t, N, D);
% signal_to_modulate2 = signal_resampled2(1:length(carrierSignal));
% DSBTC =  (1+m*signal_to_modulate2).* carrierSignal2';
% figure(7);
% plot(time, DSBTC);                             
% title('DSB-TC signal (resampled)');
% xlabel('time');
% ylabel('amplitude');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     6       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(envelope detector)
received_envelope_DSBSC = abs(hilbert(DSBSC));
received_envelope_DSBTC = abs(hilbert(DSBTC));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   7  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A, B] = rat(Fs/(Fs*10));
received_envelope_DSBSC_resampled = resample(received_envelope_DSBSC,A, B);
sound( real(double(received_envelope_DSBSC_resampled)), Fs);
pause(10)
figure(10);
plot(received_envelope_DSBSC);
title('recived DSB-SC in time domain');

received_envelope_DSBTC_resampled = resample(received_envelope_DSBTC,A, B);
sound( real(double(received_envelope_DSBTC_resampled)), Fs);
pause(10)
figure(11);
plot(received_envelope_DSBTC);
title('recived DSB-TC in time domain');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    8    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DSBTC_noise0=awgn(DSBTC,0,'measured');   %add white gausian noise
DSBTC_noise10=awgn(DSBTC,10,'measured'); %add white gausian noise
DSBTC_noise30=awgn(DSBTC,30,'measured'); %add white gausian noise

%DSB-TC with noise 0
received0=abs(hilbert(DSBTC_noise0));
received0_resampled = resample(received0,A, B);
sound( real(double(received0_resampled)), Fs);
pause(10)
figure(12);
plot(received0);
title('recived DSB-TC with noise 0 in time domain');

    
%DSB-TC with noise 10
received10=abs(hilbert(DSBTC_noise10));
received10_resampled = resample(received10,A, B);
sound( real(double(received10_resampled)), Fs);
pause(10)
figure(13);
plot(received10);
title('recived DSB-TC with noise 10 in time domain');

    
    
%DSB-TC with noise 30
received30=abs(hilbert(DSBTC_noise30));
received30_resampled = resample(received30,A, B);
sound( real(double(received30_resampled)), Fs);
pause(10)
figure(14);
plot(received30);
title('recived DSB-TC with noise 30 in time domain');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   9     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%for DSBSC Use coherent detection to receive the modulated signal with SNR=0, 10, 30 dB

%Generating noise signals
Signal_with_noise_1 = awgn(DSBSC,0,'measured'); %0 db SNR
Signal_with_noise_2 = awgn(DSBSC,10,'measured'); %10 db SNR
Signal_with_noise_3 = awgn(DSBSC,30,'measured'); %30 db SNR

%demodulating
[A, B] = rat(Fs/(Fs*10));
noisy_signal_demodulated_1 = Signal_with_noise_1.*carrierSignal';
noisy_signal_demodulated_1_resampled = resample(noisy_signal_demodulated_1, A, B);

noisy_signal_demodulated_2 = Signal_with_noise_2.*carrierSignal';
noisy_signal_demodulated_2_resampled = resample(noisy_signal_demodulated_2, A, B);

noisy_signal_demodulated_3 = Signal_with_noise_3.*carrierSignal';
noisy_signal_demodulated_3_resampled = resample(noisy_signal_demodulated_3, A, B);

%sound them 
sound( real(double(noisy_signal_demodulated_1_resampled)), Fs);
pause(10);
sound( real(double(noisy_signal_demodulated_2_resampled)), Fs);
pause(10);
sound( real(double(noisy_signal_demodulated_3_resampled)), Fs);
pause(10);

%plot them in time and frequency domain
figure(15);
plot(noisy_signal_demodulated_1);
title('recived DSB-SC with noise 0 in time domain');

figure(16);
plot(noisy_signal_demodulated_2);
title('recived DSB-SC with noise 10 in time domain');
 
figure(17);
plot(noisy_signal_demodulated_3);
title('recived DSB-SC with noise 30 in time domain');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      10       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%coherent detection with frequency error, F=100.1 KHz instead of 100 KHz

Ferror = 100100;
carrier_Signal_Error = cos(2*pi*Ferror*time);
demodulated_Frequency_Error = DSBSC.*carrier_Signal_Error';
figure(18);
plot(time, demodulated_Frequency_Error);
title('DSB-SC demodulated signal with frequency error : 100 Hz ');
xlabel('time');
%sound(real(double( resample(demodulated_Frequency_Error, A, B))),Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Repeat the coherent detection with phase error = 20

carrier_Signal_Error = cos(2*pi*Fc*time + degtorad(20));
demodulated_Phase_Error = DSBSC.*carrier_Signal_Error';
figure(19);
plot(time, demodulated_Phase_Error);
title('DSB-SC demodulated signal with frequency error : 20 degrees');
xlabel('time');
sound(real(double( resample(demodulated_Phase_Error, A, B))),Fs);





















