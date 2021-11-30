function [diag_n,theta] = array_factor(nelm,sep,phases,amplitudes)
%ARRAYFACTOR Summary of this function goes here
%   Detailed explanation goes here
theta = 0:0.001:pi;
xesp = (1:(nelm))*sep;
Fv = 0;

for i = 1:nelm
    Fv = Fv + amplitudes(i)*exp(1j*(2*pi*xesp(i)*cos(theta)+deg2rad(phases(i)))); 
end

diag = (sin(theta).^2).*Fv;
diag_n = 20*log10(abs(diag))-max(20*log10(abs(diag)));
theta = rad2deg(theta)-90;

%plot(rad2deg(theta)-90,diag_n);
%xlim([-45,45]);
%ylim([-40,0]);
end

