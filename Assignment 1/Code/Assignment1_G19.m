% Assignment (1): Relationship between solar radio flux F10.7 and sunspot number
% Team members: Ahmed Baza, Ahmmed El-Gohary, Mohamed Abbas, Ziad Baraka - Skoltech, 2020.
clc; clear all;

%% Importing and drawing data
data = importdata('data_group1.mat');
Year = data(:,1);   % Extract year
MSRF = data(:,3);   % Extract monthly solar radio flux at 10.7 cm data
MSN = data(:,4);    % Extract monthly sunspot number data
n = length(MSRF);   % Get the no. of monthly readings for each factor

date=(1947+(2/12):1/12:2016-(1/12));
figure(1)
scatter(MSN,MSRF,1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);

title('Main indicator of solar activity');
xlabel('Monthly Sunspot number');
ylabel('Montly Solar radio flux F10.7');
legend('Montly Solar radio flux F10.7');

%% Smoothing of monthly mean data by 13-month running mean. Plot results.
% Replacing the frst and last six elemetns with their means
M1 = mean(MSRF(1:6)); F1= mean(MSRF(n-5:n)); 
M2 = mean(MSN(1:6)); F2= mean(MSN(n-5:n)); 


MSRFs=MSRF;
MSNs=MSN;

% 13-month running mean
for i=7:n-6
    MSRFs(i) = (1/24)*(MSRF(i-6)+ MSRF(i+6))+(1/12)*sum(MSRF(i-5:i+5));
    MSNs(i) = (1/24)*(MSN(i-6)+ MSN(i+6))+(1/12)*sum(MSN(i-5:i+5));
end

MSRFs(1:6)= M1; MSRFs(n-5:n) = F1;
MSNs(1:6)= M2; MSNs(n-5:n) = F2;

MSNs = MSNs';

figure(2)
plot(date,MSRFs,'k',date,MSNs,'R');
title('Smoothed main indicator of solar activity');
xlabel('Date');
ylabel('Solar activity indicator');
legend('Montly Solar radio flux F10.7', 'Monthly Sunspot number');

%% Construction of multi-dimensional linear regression model
MSN0 = ones(1,826);
MSN2 = MSNs.^2;
MSN3 = MSNs.^3;

MSN_matrix = [MSN0' MSNs' MSN2' MSN3']; % Matrix of independent variables, regressors.
A = (MSN_matrix' * MSN_matrix)\ MSN_matrix' * MSRFs; % Vector of coefficients
% Reconstructing solar radio flux at 10.7 cm on the basis of sunspot number
MSRF_new = A(1) + A(2)*MSNs + A(3)*MSN2 + A(4)*MSN3;
figure(3)
plot(date,MSRFs,'k', date, MSRF_new);
title('solar radio flux model at 10.7 cm on the basis of sunspot number');
xlabel('Date');
ylabel('Solar radio flux at 10.7');
legend('Real', 'Model');

%% Determination of the variance of estimation error of solar radio flux at 10.7
Variance = (1/(n-1))*sum((MSRFs-MSRF_new').^(2));

