% This is an m-file that compares all of our denoising techniques

% Reading in our IR Spectra
% The first column is the associated wavelength for each transmittance. The
% second column is the % transmittance 
Spectra = csvread('2-4-6-trimethylphenol.CSV');

% There is alot of useless information in these files. Therefore we are
% going to trim the beginning of the data set
[A,B] = size(Spectra);
while Spectra(1,2) == 0
    Spectra(1,:) = [];
end

% Re-compute the size of Spectra and generate a random vecotor for noise
% incorporation
[A,B] = size(Spectra);
G = -1 + 2*rand(A,1);

% Plot the original Spectra for 2-4-6-trimethylphenol
figure;
plot(Spectra(1:A-1,1),Spectra(1:A-1,2),'linewidth',1.5);
title('2,4,6-trimethylphenol');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

% Incorporation of Noise
[Noisy, Noise] = addnoise(Spectra(:,2),G,35); 
Spectra2 = Spectra;
for i = 1:A
    if mod(i,2) == 0
        Spectra2(i,2) = Spectra(i,2) + Noise(i);
    end
end
Noisy = Spectra2(:,2);

% Plot the noisy Spectra for 2-4-6-trimethylphenol
figure;
plot(Spectra2(1:A-1,1),Noisy(1:A-1));
title('2,4,6-trimethylphenol with noise');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

% Re-define variables for clarity
Wavelength = Spectra(:,1);
OrigSpectra = Spectra(:,2);
NoiseSpectra = Noisy;

% Concatenate the Wavelength and Noisy Spectra
Spectra2 = [Wavelength, NoiseSpectra];

%% SVD Denoising
% Do SVD dednoising by keeping i singular values and calculate the PSNR
% with the original spectra
for i = 1:30
    
    % Perform the SVDdenoise
    SVDnoise = SVD_Denoise(Spectra2,i);
    
    % Calculate the PSNR and error in the L2 norm
    PSNRsvd(i) = PSNR(Spectra(:,2),SVDnoise(:,2));
    Errsvd(i) = sqrt(sum((SVDnoise(:,2) - Spectra(:,2)).^2));
    
    % Calculate the Pearson Corrrelation
    Meanorig = mean(Spectra(:,2)); Meanden = mean(SVDnoise(:,2));
    temp = dot(Spectra(:,2)-Meanorig,SVDnoise(:,2)-Meanden);
    Pearsonsvd(i) = temp/(norm(Spectra(:,2)-Meanorig)*norm(SVDnoise(:,2)-Meanden));
end

% Perform SVD denoising by keeping k singual values
SVDnoise1 = SVD_Denoise(Spectra2,2);
SVDnoise2 = SVD_Denoise(Spectra2,20);

% Plot the original and  denoised Spectra for SVD denoising
figure;
plot(Spectra(1:A-1,1),Spectra(1:A-1,2),'linewidth',1.5);
hold on
plot(SVDnoise1(1:A-1,1),SVDnoise1(1:A-1,2),'r');
%plot(SVDnoise2(:,1),SVDnoise2(:,2),'g');
title('2,4,6-trimethylphenol denoised by SVD');
legend('Original Spectra', 'De-noised Spectra');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
xlabel('Wavenumber (cm^-^1)');
ylabel('% Absorbtion');
%% MAD Denoising

% We need to chop the original Spectra vector for comparison
MADorigSpectra(:,1) = ChopVector(Spectra(:,1),6);
MADorigSpectra(:,2) = ChopVector(Spectra(:,2),6);

% Perform MAD denoising with the Daub(h) filter with its iterations and
% calculate the PSNR
for its = 1: 6
    for h = 2:2:6
        
        % Perfrom the MAD_Denoise Algorithm
        MADnoise = MAD_Denoise(Spectra2,Daub(h),its);
        
        % Calculate the PSNR and error with the L2 norm
        PSNRmad(its,h) = PSNR(MADorigSpectra(:,2),MADnoise(:,2));
        errmad(its,h) = sqrt(sum((MADorigSpectra(:,2)-MADnoise(:,2)).^2));
        
        % Calculate the Pearson Corrrelation
        Meanorig = mean(MADorigSpectra(:,2)); Meanden = mean(MADnoise(:,2));
        temp = dot(MADorigSpectra(:,2)-Meanorig,MADnoise(:,2)-Meanden);
        Pearsonmad(its,h) = temp/(norm(MADorigSpectra(:,2)-Meanorig)*norm(MADnoise(:,2)-Meanden));
    end
end

MADnoise = MAD_Denoise(Spectra2,Daub(4),4);
% Plot the denoised Spectra for MAD denoising
figure;
plot(MADorigSpectra(:,1),MADorigSpectra(:,2),'linewidth',1.5);
hold on
plot(MADnoise(:,1),MADnoise(:,2),'r');
title('2,4,6-trimethylphenol denoised by MAD');
legend('Original Spectra', 'De-noised Spectra');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
xlabel('Wavenumber (cm^-^1)');
ylabel('% Absorbtion');
%% SURE Denoising

% We need to chop the original Spectra vector for comparison
SUREorigSpectra(:,1) = ChopVector(Spectra(:,1),6);
SUREorigSpectra(:,2) = ChopVector(Spectra(:,2),6);

% Perform SURE denoising with the Daub(h) filter with its iterations and
% calculate the PSNR and the error
for its = 1: 6
    for h = 2:2:6
        
        % Perfrom the SURE_Denoise Algorithm
        SUREnoise = SURE_Denoise(Spectra2,Daub(h),its);
        
        % Calculate the PSNR and error with the L2 norm
        PSNRsure(its,h) = PSNR(SUREorigSpectra(:,2),SUREnoise(:,2));
        errsure(its,h) = sqrt(sum((SUREorigSpectra(:,2) - SUREnoise(:,2)).^2));
    
        % Calculate the Pearson Corrrelation
        Meanorig = mean(SUREorigSpectra(:,2)); Meanden = mean(SUREnoise(:,2));
        temp = dot(SUREorigSpectra(:,2)-Meanorig,SUREnoise(:,2)-Meanden);
        Pearsonsure(its,h) = temp/(norm(SUREorigSpectra(:,2)-Meanorig)*norm(SUREnoise(:,2)-Meanden));
    end
end

SUREnoise = SURE_Denoise(Spectra2,Daub(6),4);
% Plot the denoised Spectra for SURE denoising
figure;
plot(SUREorigSpectra(:,1),SUREorigSpectra(:,2),'linewidth',1.5);
hold on
plot(SUREnoise(:,1),SUREnoise(:,2),'r');
title('2,4,6-trimethylphenol denoised by SURE');
legend('Original Spectra', 'De-noised Spectra');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
xlabel('Wavenumber (cm^-^1)');
ylabel('% Absorbtion');