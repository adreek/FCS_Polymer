function [ G ] = fcspddeb( x, tau )
% Primary Equation: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.73.041919
% Correlation Function which incorporated intermolecular dynamics of linearized DNA in FCS on Matlab.

t=tau
S=5
F = 0.00634639161199272
tau_f = 0.005
% G = ((1-x(4) + x(4)*exp(-tau./0.005))./(1-x(4))).*(1/x(1)*((1+tau/x(3)).^-1).*(1+tau/(8^2*x(3))).^(-1/2)) + 1
f=8

% g_bar = 1e-7 %scaling factor  x_1
% L = 4500 %flexibility coeff x_2
% f = 7 %related to aspect ratio x_3
% r = 3 %waist size of beam (nm) x_4
% D = 1e-2%diff coeff of int dyn (cm^2/s) x_5
% gamma = 0.05 %tunable decay rate x_6

% n is N
% d is diffusion time in ms
% g is gamma

x1_n = 1
x1_d = 50.4887 %Input Diffusiob Time (micro)
x1_g = 0.0849

x2_n = 1
x2_d = 95.1184
x2_g = 0.0603

G = ((1-F + F*exp(-tau./0.005))./(1-F)).*(x(1)*((((1+tau/x1_d).^-1).*(1+tau/(f.^2*x1_d)).^(-1/2)).*((((2*x1_n*1E6*f)*(sqrt((3*2.5e-7.^2)/(2*(f.^2-1)))))*(1 + 1./x1_d + 1.65*((x1_g*t).^(3/4))).^(-1/2))).*atan((sqrt((2*(f.^2-1)*1E6^2)/(3*2.5e-7.^2))).*(1 + 1./x1_d + 1.65*(x1_g*t).^(3/4)).*(f.^2+(2*1E6.^2/(3*2.5e-7.^2))+ 1./x1_d + 1.65*(x1_g*t).^(3/4)).^(-1/2)))+x(2)*(((1+tau/x2_d).^-1).*(1+tau/(f.^2*x2_d)).^(-1/2)).*((((2*x2_n*1E6*f)*(sqrt((3*2.5e-7.^2)/(2*(f.^2-1)))))*(1 + 1./x2_d + 1.65*((x2_g*t).^(3/4))).^(-1/2))).*atan((sqrt((2*(f.^2-1)*1E6^2)/(3*2.5e-7.^2))).*(1 + 1./x2_d + 1.65*(x2_g*t).^(3/4)).*(f.^2+(2*1E6.^2/(3*2.5e-7.^2))+ 1./x2_d + 1.65*(x2_g*t).^(3/4)).^(-1/2)))
% 
end
