%Algoritmo de filtros low-pass FIR e IIR
%Autor: André Heidemann Iarozinski
%Data: 28/11/2015

fs=44100; % definindo freq de amostragem
recc = audiorecorder(44100, 16, 1);
disp('inicio da gravacao')
recordblocking(recc, 5); % duração da gravação: 5 seg
disp('fim da gravacao');
u = getaudiodata(recc);
S = fft(u);
L = length(S); % normalizando S
P2 = abs(S/L);
P = P2(1:L/2+1);
P(2:end-1) = 2*P(2:end-1);
f = fs*(0:(L/2))/L;
plot(f,P) % plotando a FFT de 0Hz até fs/2 da amostra de áudio


%% Filtro passa baixas FIR
fa = 44100; % frequencia de amost
s1 = u; %% sinal gravado
fp = 7000; % frequência de passagem
fs = 7500; % frequência de corte
% normalização das frequências
wp = (fp/(fa/2))*pi
ws = (fs/(fa/2))*pi
bt = ws - wp; %banda de transição
M = ceil((6.6*pi/bt)) + 1; % M de acordo com a tabela das janelas
wc = (ws + wp)/2; %frequência de corte intermediária
alfa = (M-1)/2; %% filtro passa baixas ideal
n=0:M-1;
m = n- alfa + eps;
hd = sin(wc*m) ./ (pi*m) ; % resposta impulsiva do fpb ideal
jan = hamming(M)'; %calcula a janela de hamming
h = hd.*jan; % multiplicação entre os vetores
sinal_filtrado = conv(h,s1); %convolução entre os sinais
sound(sinal_filtrado,fa);
S = fft(sinal_filtrado);
L = length(S); % normalizando S
P2 = abs(S/L);
P = P2(1:L/2+1);
P(2:end-1) = 2*P(2:end-1);

f = fa*(0:(L/2))/L;
plot(f,P) % plotando a FFT do sinal filtrado de 0Hz até fs/2


%% filtro passa-baixa IIR butterworth
rp=2;
as=30;
wp=0.25*pi; %%% freq de borda da banda de passagem
ws=0.35*pi; %%% freq de borda da banda de rejeição
T=1;
wap = wp/T;
was = ws/T;
%% prototipo do filtro anlogico 
N = ceil( log10 ( (10^(rp/10) -1) / (10^(as/10) -1) ) / (2*log10(wp/ws)));
wc = (wp/(((10^(rp/10)-1))^(1/(2*N))));
[z,p,k] = buttap(N); %% retorna os zeros, os polos e o ganho
num = real(poly(z)); % numerador da func de transf Ha(s)
num = num*(wc^N)*k;
den= real(poly(p*wc)); %%% denomidador da funcao de transf
[numd,dend] = impinvar(num,den,T); %% discretrizando
sys =tf(numd,dend);
x=filter(numd,dend,u); % filtragem do sinal
sound(x,44100);
w=0:pi/100:pi;
H = freqz(numd,dend,w); %% Ha(s)
Hma = abs(H);
Hfase = angle(H);
Hmdb = 20*log10((Hma+eps)/(max(Hma)));
plot(w/pi, Hmdb); % filtro em dB

