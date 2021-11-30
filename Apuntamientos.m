%% Common data
clc;
clear;
close all;

nelm = 8;
separation = 0.55; %Lambdas
%% Aims
aims = [5,10,15,30];
phases = zeros(length(aims),nelm);
for i = 1:length(aims)
    prog_phase(i) = generate_aim(aims(i),separation);
    phases(i,:) = wrapTo360(prog_phase(i).*(0:(nelm-1)));
    [phases_d(i,:)] = discretize_phase(phases(i,:),6);  
end
writePortConfigCST("Array 8x8 - elemento stub - polV - fabr - Copy.cst",...
    ones(1,nelm),phases(2,:),["beam10"]);
figure(1);
for i = 1:length(aims)
    [diagram,theta] = array_factor(nelm,separation,phases_d(i,:),ones(1,nelm));
    plot(theta,diagram);
    if(i == 1)
        hold on;
    end
end
hold off;
xlim([-45,45]);
ylim([-40,0]);
legend(["5","10","15","30"]);
%% Amplitudes
amps = zeros(5,nelm);
amps_d = zeros(5,nelm);
titles = ["Coseno sobre pedestal","Cuadrada","Triangular","Chebyshev","Taylor"];
%Coseno_SobrePedestal
N = 10;
H = 0.4;
pos = linspace(-(nelm-1)/2,(nelm-1)/2,nelm);
amps(1,:) = 1+H*cos(pi*pos/(N-1)).^2;
amps(1,:) = 20*log10(amps(1,:))-max(20*log10(amps(1,:)));
[amps_d(1,:)] = discretize_amplitude(amps(1,:),6,0.5);
%Squared_distribution
delta = 0.5;
amps(2,:) = 1-(1-delta)*(pos/3.5).^2;
amps(2,:) = 20*log10(amps(2,:))-max(20*log10(amps(2,:)));
[amps_d(2,:)] = discretize_amplitude(amps(2,:),6,0.5);
%Triangular
amps(3,1:4) = 1+2*pos(1:4)/8;
amps(3,5:8) = flip(amps(3,1:4));
amps(3,:) = 20*log10(amps(3,:))-max(20*log10(amps(3,:)));
[amps_d(3,:)] = discretize_amplitude(amps(3,:),6,0.5);
%Cheb
amps(4,:) = chebwin(nelm,20);
amps(4,:) = 20*log10(amps(4,:))-max(20*log10(amps(4,:)));
[amps_d(4,:)] = discretize_amplitude(amps(4,:),6,0.5);
%Taylor
amps(5,:) = taylorwin(nelm,5,-20);
amps(5,:) = 20*log10(amps(5,:))-max(20*log10(amps(5,:)));
[amps_d(5,:)] = discretize_amplitude(amps(5,:),6,0.5);

figure(2);
for i = 1:5
    subplot(3,2,i)
    plot(pos,amps(i,:),"*");
    hold on;
    plot(pos,amps_d(i,:),"*");
    title(titles(i));
    hold off;
end

figure(3);
for i = 1:5
    subplot(3,2,i)
    [diagram,theta] = array_factor(nelm,separation,zeros(1,nelm),10.^(amps(i,:)/20));
    plot(theta,diagram);
    hold on;
    [diagram,theta] = array_factor(nelm,separation,zeros(1,nelm),10.^(amps_d(i,:)/20));
    plot(theta,diagram);
    xlim([-45,45]);
    ylim([-40,0]);
    title(titles(i));
    hold off;
end

%% Amplitude distributions

%% Amplitude and aim

function [prog_phase] = generate_aim(direction,separation)
    prog_phase = rad2deg(2*pi*separation*sind(direction));
end

function [amp_val,amp_chip] = discretize_amplitude(amplitudes,n_bits,resolution)
    bits = 0:(2^n_bits-1);
    scale = bits*(-resolution);
    amp_chip = zeros(size(amplitudes));
    amp_val = zeros(size(amplitudes));
    for i = 1:length(amplitudes)
        [~,amp_chip(i)] = min(abs(scale-amplitudes(i)));
        amp_val(i) = scale(amp_chip(i));
    end
end


function [pha_val,pha_chip] = discretize_phase(phases,n_bits)
    scale = 0:(360/(2^n_bits)):(360-360/(2^n_bits));
    pha_chip = zeros(size(phases));
    pha_val = zeros(size(phases));
    for i = 1:length(phases)
        [~,pha_chip(i)] = min(abs(scale-phases(i)));
        pha_val(i) = scale(pha_chip(i));
    end
end

