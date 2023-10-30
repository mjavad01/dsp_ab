%%2.1.a
x = [1, 2, 3];
h = [4, 3, 2, 1];
myconv(x,h)
conv(x,h)


%%2.1.b

n = 0:200;
k = 50;
f = 1/k;
h = zeros(1, length(n));
h(1:10) = 0.1 * 1;
x = square(2 * pi * f * n);
y = myconv(h, x);

figure(1);
subplot(3,1,1);
stem(n, h,'g');
xlabel('n') ;
ylabel('') ;
title('H[n]') ;
grid on ;

subplot(3,1,2);
stem(n, x,'r');
xlabel('n') ;
ylabel('') ;
title('X[n]') ;
grid on ;

subplot(3,1,3);
stem(0:length(y) - 1, y,'m');
xlabel('n') ;
ylabel('') ;
title('Y[n]') ;
grid on ;


%2.1.c

n = 0:14;
h = 0.25 * 0.75.^n;
y = myconv(h, x);

figure(2);
subplot(3,1,1);
stem(n, h,'g');
xlabel('n') ;
ylabel('') ;
title('H[n]') ;
grid on ;

subplot(3,1,2);
stem(0:length(x) - 1, x,'r');
xlabel('n') ;
ylabel('') ;
title('X[n]') ;
grid on ;

subplot(3,1,3);
stem(0:length(y) - 1, y,'c');
xlabel('n') ;
ylabel('') ;
title('Y[n]') ;
grid on ;


%2.1.d

hn = 0.2*[1,-5,10,-10,5,-1];
y = myconv(hn, x);

figure(3);
subplot(3,1,1);
stem(0:length(hn) - 1, hn,'g');
xlabel('n') ;
ylabel('') ;
title('H[n]') ;
grid on ;

subplot(3,1,2);
stem(0:length(x) - 1, x,'r');
xlabel('n') ;
ylabel('') ;
title('X[n]') ;
grid on ;

subplot(3,1,3);
stem(0:length(y) - 1, y,'m');
xlabel('n') ;
ylabel('') ;
title('Y[n]') ;
grid on ;


%2.2.a

M  = 100;
n  = 0 : M - 1;
w1 = 0.05 * pi;
w2 = 0.20 * pi;
w3 = 0.35 * pi;
wa = 0.15 * pi;
wb = 0.25 * pi;


s  = sin(w2 * n);
v  = sin(w1 * n) + sin (w3 * n);
x  = s + v;

w  = 0.54 - 0.46 * sin(2 * pi * n / M);
h  = w .* ( wb * sinc(wb * (n - M/2)/pi)/pi - wa * sinc(wa * (n - M/2)/pi)/pi );

figure(4);
subplot(2,1,1);
stem(n, x,'g');
xlabel('n') ;
ylabel('') ;
title('X[n]') ;
grid on ;

subplot(2,1,2);
stem(n, s,'m');
xlabel('n') ;
ylabel('') ;
title('S[n]') ;
grid on ;


%2.2.b

y = filter(h,1,x);

figure(5);
subplot(2,1,1);
stem(n, s,'g');
xlabel('n') ;
ylabel('') ;
title('S[n]') ;
grid on ;

subplot(2,1,2);
stem(n, y,'b');
xlabel('n') ;
ylabel('') ;
title('Y[n]') ;
grid on ;


%2.2.c

load('Filter1.mat');
y = filter(Num,1,x);
 
figure(6);
subplot(2,1,1);
stem(n, s,'g');
xlabel('n') ;
ylabel('') ;
title('S[n]') ;
grid on ;

subplot(2,1,2);
stem(n, y,'b');
xlabel('n') ;
ylabel('') ;
title('Y[n]') ;
grid on ;


%2.3.a

[x, Fs] = audioread('Audio01.wav');
T = 10;
t = 0: 1/Fs:T-1/Fs;
L_x       = length(x);
f_x       = (Fs/L_x) * (-L_x/2:L_x/2-1);
fft_x     = fftshift(fft(x))/L_x;

figure(7);
subplot(2, 1, 1);
plot(t, x, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time(Sec)') ;
ylabel('Amplitude') ;
title('Signal in Time Domain');

subplot(2,1,2) ;
plot(f_x, abs(fft_x), 'r', 'LineWidth',1.5) ;
grid on;
xlabel('freq (Hz)');
ylabel('|Signal (f)|');
title('Signal in Frequency Domain');


%2.3.b

load('Filter2.mat');
y_0 = filter(Num2,1,x);

% L_y_0       = length(y_0);
% f_y_0       = (Fs/L_y_0) * (-L_y_0/2:L_y_0/2-1);
% fft_y_0     = fftshift(fft(y_0))/L_y_0;
% 
% figure(99);
% subplot(2, 1, 1);
% plot(t, y_0, 'g', 'LineWidth', 1.5);
% grid on;
% xlabel('Time(Sec)') ;
% ylabel('Amplitude') ;
% title('Signal in Time Domain');
% 
% subplot(2,1,2) ;
% plot(f_y_0, abs(fft_y_0), 'r', 'LineWidth',1.5) ;
% grid on;
% xlabel('freq (Hz)');
% ylabel('|Signal (f)|');
% title('Signal in Frequency Domain');


%2.3.c

f0 = 10000;
n   = 1:length(x);
s_n = 2*cos( 2 * pi * f0 * n / Fs );

y_1 = y_0 .* s_n';
y_2 = filter(Num2,1,y_1);

% soundsc(y_2, Fs);


%2.3.d

y_3 = filter(Num2,1,y_2);
y_4 = y_3 .* s_n';
Sig_filtered = filter(Num,1,y_4);
soundsc(Sig_filtered, Fs);


%2.4.a

f0    = 100;
fs    = 500;
t_min = 0;
t_max = 2;
t     = t_min : 1/fs : t_max - 1/fs;
s_1   = sin(2 * pi * f0 * t);
f0    = 200;
f1    = 400;
s_2   = chirp(t,f0,t_max,f1);
s_3   = zeros(1,length(s_1));
s_3(250) = 50;
s     = s_1 + s_2 + s_3;

% figure();
% subplot(4,1,1);
% plot(t, s_1);
% subplot(4,1,2);
% plot(t, s_2);
% subplot(4,1,3);
% plot(t, s_3);
% subplot(4,1,4);
% plot(t, s);


L_x       = length(s);
f_x       = (fs/L_x) * (-L_x/2:L_x/2-1);
fft_x     = fftshift(fft(s))/L_x;

figure(8);
plot(f_x, abs(fft_x), 'r', 'LineWidth',1.5) ;
grid on;
xlabel('freq (Hz)');
ylabel('|Signal (f)|');
title('Signal in Frequency Domain');


%2.4.b

L = [256, 512];
w = hamming(L(1)/2);
figure(9);
spectrogram (s, w);

w = hamming(L(2)/2);
figure(10);
spectrogram (s, w);


%2.5

loc = linspace(0,1,2^10);
[x,x_noisy] = wnoise('doppler',10,7);

figure(11);
subplot(2,1,1);
plot(loc,x, 'g', 'LineWidth',1.5);
grid on;
title('Clean Doppler');

subplot(2,1,2);
plot(loc,x_noisy, 'c', 'LineWidth',1.5);
grid on;
title('Noisy Doppler');

[cA1,cD1] = dwt(x_noisy,'db1');
loc = linspace(0,1,length(cA1));

figure(12);
subplot(3,1,1);
plot(loc,cA1, 'r', 'LineWidth',1.5);
grid on;
title('Level 1');

[cA2,cD2] = dwt(cA1,'db1');
loc = linspace(0,1,length(cA2));
subplot(3,1,2);
plot(loc,cA2, 'm', 'LineWidth',1.5);
grid on;
title('Level 2');

[cA3,cD3] = dwt(cA2,'db1');
loc = linspace(0,1,length(cA3));
subplot(3,1,3);
plot(loc,cA3, 'g', 'LineWidth',1.5);
grid on;
title('Level 3');


function y = myconv(h,x)
M = length(h);
L = length(x);
N = M + L - 1;
X = [x, zeros(1,M)];
H = [h, zeros(1,L)];
y = zeros(1, N);
if(M > L)
    temp = X;
    X    = H;
    H    = temp;
end
for ii = 1 : N
    for jj = 1 : min(L, M)
        if (ii - jj + 1 > 0)
            y(ii) = y(ii) + H(jj) * X(ii - jj + 1);
        end
    end
end
end


