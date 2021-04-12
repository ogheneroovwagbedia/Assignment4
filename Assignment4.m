%% oghenero ovwagbedia
%% 101040228

% Parameters

R1 = 1;
Ca = 0.25;
R2 = 2;
L = 0.2;
R3 = 10;
alpha = 100;
R4 = 0.1; 
R0 = 1000;
% making matrices
G0 = 1/R0;
G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;

G=zeros(8);
C=zeros(8);

%G matrix
G(1,:)=[G1 -G1 0 0 0 0 0 G1];
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/R0) 0 0 0]; 
G(6,:)=[0 G1 -G1 0 0 0 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
G(8,:)=[1 0 0 0 0 0 0 0]; 
%C matrix
C(1,:)=[Ca -Ca 0 0 0 0 0 0]; 
C(2,:)=[-Ca Ca 0 0 0 0 0 0];
C(3,:)=[0 0 0 0 0 0 0 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
C(6,:)=[0 0 0 0 0 -L 0 0]; 
C(7,:)=[0 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0]; 

V1 = [];
V2 = [];

for Vin=-10:1:10
    F=[0; 0; 0; 0; 0; 0; 0; Vin];
    V=G\F;
    V1 = [V1 V(1)];
    V2 = [V2 V(5)];
end

figure(1)
hold on;
title('DC Sweep');
xlabel('Vin sweep from -10 to 10 (V)');
ylabel('Voltage (V)');
plot(-10:1:10, V1);
plot(-10:1:10, V2);
hold off;
legend('V3', 'VO');
count=1000;
om=zeros(2,count);
om(1,:)=linspace(0,500,count);
Vin=1;
for i=1:count
    og=om(1,i);
    F=[0; 0; 0; 0; 0; 0; 0; Vin];
    V=(G+1j*og*C)\F;
    om(2,i)=V(5);
end

figure(2)
plot(om(1,:),real(om(2,:)));
title('AC plot - VO as a function of Omega');
ylabel('VO (V)');
xlabel('radians/s');
figure(3)
V2 = [];
w = 3.14;
Vin = 1;
F=[0; 0; 0; 0; 0; 0; 0; Vin];

for w=1:1:10
    ep = (G+2*w^2*1j*C)\F;
    V2 = [V2 20*log10(abs(ep(5)/F(8)))];
end

semilogx(1:1:10, V2);
hold on;
title('AC Sweep');
xlabel('Radians/sec');
ylabel('dB');
grid on;
% figure(4);
 cmd =  Ca + 0.05.*randn(5000,1) ;
for i = 1:100
    F = [0 0 0 0 0 1 0 0];
    CMND =   [ -cmd(i)  cmd(i) 0  0  0  0  0  0;
                cmd(i) -cmd(i) 0  0  0  0  0  0;
                0  0 0  0  0  0  0  0; 
                0  0 0  0  0  0  0  0;
                0  0 0  0  0  0  0  0; 
                0  0 0  0  0  0  0  0; 
                0  0 0  0  0  0  L  0; 
                0  0 0  0  0  0  0  0];
    V = (G+(pi*CMND))\F';
    Vh(i) = V(5);
end
figure(4)
hist(Vh)
title('Histogram of Gain')
%% part 2
step = 1000;
Cn = 0;
t = 1
I = zeros(1,step);
dstep = t/step;

Vin = zeros(1,step);
Vin(0.03*step:step) = 1;
F = zeros(8,1,step);

for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

figure(5);
Vout = k(5,:,:);
Vout = Vout(1,:);
hold on;
plot(linspace(0,t,step), Vout);
plot(linspace(0,t,step), Vin);
title('Step input voltages');
legend('Vout', 'Vin');
xlabel('Time (s)');
ylabel('V (V)');

Vin = sin(linspace(0,1,step)*2*pi*1/0.03);
F = zeros(8,1,step);

for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

Vout = k(5,:,:);
Vout = Vout(1,:);
figure(6);
hold on;
plot(linspace(0,t,step), Vout);
plot(linspace(0,t,step), Vin);
title('Voltage with sine input');
legend('Vout', 'Vin');
xlabel('Time (s)');
ylabel('Vo (V)');

Vin = gaussmf(linspace(0,1,step),[0.03 0.06]);

F = zeros(8,1,step);
for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

Vout = k(5,:,:);
Vout = Vout(1,:);
figure(7);
hold on;
plot(linspace(0,t,step), Vout);
plot(linspace(0,t,step), Vin);
hold off;
xlabel('Time (s)');
ylabel('Vout (V)');
title('Voltage for Gaussian Function');
legend('Vout', 'Vin');

Vin = sin(linspace(0,1,step)*2*pi*1/0.03);
F = zeros(8,1,step);

for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

F = abs(fftshift(fft(Vout)));
figure(8);
hold on;
plot(((1:length(F))/step)-0.5,20*log10(F));

F = abs(fftshift(fft(Vin)));
plot(((1:length(F))/step)-0.5,20*log10(F));

xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Sine Function Frequency Response');

Vin = gaussmf(linspace(0,1,step),[0.03 0.06]);

F = zeros(8,1,step);
for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

F = abs(fftshift(fft(Vout)));

figure(9);
hold on;
plot(((1:length(F))/step)-0.5,20*log10(F));
F = abs(fftshift(fft(Vin)));
plot(((1:length(F))/step)-0.5,20*log10(F));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Gaussian Function Frequency Response');

%% part 3
% it can be seen that increasing the value of the capacitance also increases the effect of noise on the circuits
% frequency response.
% updated C matrix 

C(1,:)=[Ca -Ca 0 0 0 0 0 0]; 
C(2,:)=[-Ca Ca 0 0 0 0 0 0];
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
C(6,:)=[0 0 0 0 0 -L 0 0]; 
C(7,:)=[0 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0]; 


Vin = gaussmf(linspace(0,1,step),[0.03 0.06]);
I = 0.001*rand(step,1);

F = zeros(8,1,step);
for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

Vout = k(5,:,:);
Vout = Vout(1,:);
figure(10);
hold on;
plot(linspace(0,t,step), Vout);
plot(linspace(0,t,step), Vin);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Voltages for Gassian plot with noise');
legend('Vout', 'Vin');

FN = abs(fftshift(fft(Vout)));
figure(11);
hold on;
plot(((1:length(FN))/step)-0.5,20*log10(FN));
FN = abs(fftshift(fft(Vin)));
plot(((1:length(FN))/step)-0.5,20*log10(FN));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Gaussian pluse with noise');

figure(12);
hold on;
FN = abs(fftshift(fft(Vout)));
plot(((1:length(FN))/step)-0.5,20*log10(FN));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
title('Voltages with noise with different Cn');

Cn = 0.0001;

G(1,:)=[1 -1 0 0 0 0 0 1];
C(1,:)=[Ca -Ca 0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0]; 
C(2,:)=[-Ca Ca 0 0 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/R0) 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
G(6,:)=[0 1 -1 0 0 0 0 0]; 
C(6,:)=[0 0 0 0 0 -L 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
C(7,:)=[0 0 0 0 0 0 0 0]; 
G(8,:)=[1 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0];

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

Cn = 0.01;

G(1,:)=[1 -1 0 0 0 0 0 1];
C(1,:)=[Ca -Ca 0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0]; 
C(2,:)=[-Ca Ca 0 0 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/R0) 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
G(6,:)=[0 1 -1 0 0 0 0 0]; 
C(6,:)=[0 0 0 0 0 -L 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
C(7,:)=[0 0 0 0 0 0 0 0]; 
G(8,:)=[1 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0];

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

Vout = k(5,:,:);
Vout = Vout(1,:);
FN = abs(fftshift(fft(Vout)));
plot(((1:length(FN))/step)-0.5,20*log10(FN));

Cn = 0.00000001;

G(1,:)=[1 -1 0 0 0 0 0 1];
C(1,:)=[Ca -Ca 0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0]; 
C(2,:)=[-Ca Ca 0 0 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/R0) 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
G(6,:)=[0 1 -1 0 0 0 0 0]; 
C(6,:)=[0 0 0 0 0 -L 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
C(7,:)=[0 0 0 0 0 0 0 0]; 
G(8,:)=[1 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0];

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

Vout = k(5,:,:);
Vout = Vout(1,:);
FN = abs(fftshift(fft(Vout)));
plot(((1:length(FN))/step)-0.5,20*log10(FN));
legend('Cn = 0.1', 'Cn = 0.01', 'Cn = 0.0000000001');



t = 1;

dstep = t/step;
Vin = gaussmf(linspace(0,1,step),[0.03 0.06]);
I = 0.001*rand(step,1);

F = zeros(8,1,step);
for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

Vout = k(5,:,:);
Vout = Vout(1,:);

figure(13);
hold on;
subplot(3,1,1);
plot(linspace(0,t,step), Vout);
xlabel('t (s)');
ylabel('Vout (V)');
title('Gaussian Function Transient Response with Various Time Steps - 100');

step = 1000;
t = 1;
tau= 0.03;
mean = 0.06;
dstep = t/step;
Vin = gaussmf(linspace(0,1,step),[tau mean]);
I = 0.001*rand(step,1);


F = zeros(8,1,step);
for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end
Vout = k(5,:,:);
Vout = Vout(1,:);
subplot(3,1,2);
plot(linspace(0,t,step), Vout);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Gaussian Function Transient Response with Various Time Steps - 1000');

step = 100000;
t = 1;
dstep = t/step;
Vin = gaussmf(linspace(0,1,step),[tau mean]);
I = 0.001*rand(step,1);

F = zeros(8,1,step);
for i=1:step
    F(3,1,i) = -I(i);
    F(8,1,i) = Vin(i);
end

k = zeros(8,1, step);

for i=2:step
    ep = C/dstep + G;
    k(:,:,i) = ep\(C*k(:,:,i-1)/dstep + F(:,:,i));
end

Vout = k(5,:,:);
Vout = Vout(1,:);

subplot(3,1,3);
plot(linspace(0,t,step), Vout);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Gaussian Function Transient Response with Various Time Steps - 100000');

%% question 6 NON -linearity comments
%In order to implement this a new column matrix has to be introduced into
%the equation(Gx + f(x) = b). And to solve the system the jacobian of that
%matrix will need to be solved.
%For the simulation of the new equation of the circuit, the inclusion of
%that equation would increase the size of the matrix.




