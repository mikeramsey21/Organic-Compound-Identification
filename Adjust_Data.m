% This is a Matlab Script to do preliminary work with my IR Spectroscopy
% data. I obatained data for 9 different samples and I attend to perform
% denoising and create a data bank for comparison of Spectra. This
% particular file will be used to adjust the IR data so that the baseline
% is at 100% absobance

% Reading in our IR Spectra
% The first column is the associated wavelength for each transmittance. The
% next columns are the % transmittance for each compound
Spectra = csvread('IR_Spectra_Combined.csv',1);
[A,B] = size(Spectra);

for i = 2:B
    if Spectra(A-1,i) < 100
        scale = 100 - Spectra(A-1,i);
        Spectra(261:A,i) = Spectra(261:A,i) + scale;
    end
    for j = 261:A
        if Spectra(j,i) > 100
            Spectra(j,i) = 100;
        end
    end
end

csvwrite('IR_Spectra_Adjusted.csv',Spectra);