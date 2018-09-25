% This is an m-file to perfrom functional group detection

% Reading in our IR Spectra
% The first column is the associated wavelength for each transmittance. The
% second column is the % transmittance 
Spectra = csvread('4-nitrobenzaldehyde.CSV');

% There is alot of useless information in these files. Therefore we are
% going to trim the beginning of the data set
[A,B] = size(Spectra);
while Spectra(1,2) == 0
    Spectra(1,:) = [];
end
[A,B] = size(Spectra);
Spectra(A,:) = [];

% Re-compute the size of Spectra and generate a random vecotor for noise
% incorporation
[A,B] = size(Spectra);
% 
% % Perform naive baseline correction
% scale = 100 - Spectra(A,2);
% for i = 1:A
%     if Spectra(i,2) < 100
%         Spectra(i,2) = Spectra(i,2) + scale;
%     end
%     if Spectra(i,2) >= 100
%         Spectra(i,2) = 100;
%     end
% end

% Perform n iterations of Signal Smoothing
% The divisors of 7208
% 1, 2, 4, 8, 17, 34, 53, 68, 106, 136, 212, 424, 901, 1802, 3604, 7208
n1 = 5; n3 = 20; n = n1 + n3;
its1 = Smooth(Spectra,150,n1);
its3 = Smooth(its1,34,n3);

% Plot the result
figure;
plot(its3(:,1),its3(:,2),'linewidth',1.5);
hold on
plot(its3(:,2*n+1),its3(:,2*n+2),'linewidth',1.5);
title('4-nitrobenzaldehyde');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
xlabel('Wavenumber (cm^-^1)');
ylabel('% Absorbtion');
hold on

% Perform peak detection
Spectra2 = [its3(:,2*n+1),its3(:,2*n+2)];
Npeaks = peaks(Spectra2,300);

% Plot the peaks
scatter(Npeaks(:,1),Npeaks(:,2),'k','linewidth',2);
legend('Original Spectra','Smoothed Spectra','Peaks');

% Find out the function groups that are present
Identify(Npeaks)
