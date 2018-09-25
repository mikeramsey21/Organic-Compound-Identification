% This is a Matlab Script to do preliminary work with my IR Spectroscopy
% data. I obatained data for 9 different samples and I attend to perform
% denoising and create a data bank for comparison of Spectra. The latter
% will be used for identification of compounds. 

% Reading in our IR Spectra
% The first column is the associated wavelength for each transmittance. The
% next columns are the % transmittance for each compound
Spectra = csvread('IR_Spectra_Combined.csv',1);

% There is alot of useless information in these files. Therefore we are
% going to trim the beginning of the data set
[A,B] = size(Spectra);
while Spectra(1,2:B) == zeros(1,B-1)
    Spectra(1,:) = [];
end

% The length of our vectors are odd.  We can use the DiscreteWavelets
% Toolbox command |ChopVector| to remove an elements from the end of
% Spectra to create vectors with a length divisible by 64. We need this to
% utilize multiple iterations of Daubechies D4 filter
for i = 1:B
    Spectra2(:,i) = ChopVector(Spectra(:,i),6);
end
[A,B] = size(Spectra2);

%% Plot the Spectra for 2-4-6-trimethylphenol
figure;
plot(Spectra2(:,1),Spectra2(:,2),'linewidth',1.5);
title('2,4,6-trimethylphenol');
grid on; grid minor;
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
xlabel('Wavenumber (cm^-^1)');
ylabel('Percent absorbtion');

%%
% Plot one iteration iteration
its = 1;
h = Daub(4)
Spectra4 = WT1D(Spectra2(:,2),h,its);

% Plot it
figure;
WaveletVectorPlot(Spectra4,1);

%% We will start with three interations of the D4 filter on each column vector
its = 3;
h = Daub(4);
for i = 1:B
Spectra3(:,i) = WT1D(Spectra2(:,i),h,its);
end

% Plot 3 iterations of the wavelet transform for 2-4-6-trimethylphenol
figure;
WaveletVectorPlot(Spectra3(:,2),3);

% Try four interations of the D4 filter on each column vector
its = 4;
h = Daub(4);
for i = 1:B
Spectra4(:,i) = WT1D(Spectra2(:,i),h,its);
end

% Plot 4 iterations of the wavelet transform for 2-4-6-trimethylphenol
figure;
WaveletVectorPlot(Spectra4(:,2),4);

% Compute the cumulative energy of the original spectrum and the transforms
ce1 = CE(Spectra2(:,2));
ce2 = CE(Spectra3(:,2));
ce3 = CE(Spectra4(:,2));

% Plot the various cumulative energies
figure;
plot(ce1);
hold on
plot(ce2,'r');
plot(ce3,'g');
title('Cumulative energy of the WT');
legend('Original Spectra','3 iterations','4 iterations'); 
hold off;

%  We will retain 99.99% of the energy of the 3 iteration transform.  We can use the
%  command |nCE| to determine how many elements of |Spectra3(:,1)| contribute to this
%  energy level.
k1 = nCE(ce2, 0.9999);
disp(sprintf('99.99%% of the energy is stored in the largest %i (in modulus) elements.', k1));
disp(sprintf('We convert %i elements in B to 0.',length(Spectra3(:,2))-k1));
newSpectra3(:,1) = Comp(Spectra3(:,2),k1);

% Calculate the inverse transform
its = 3;
compressedSpectra = IWT1D(newSpectra3,h,its);

% Re-plot the spectra for 2-4-6-trimethylphenol after compression
figure;
plot(Spectra2(:,1),compressedSpectra);
title('Compressed 2,4,6-trimethylphenol');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

