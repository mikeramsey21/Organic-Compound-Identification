function [denSpectra] = SURE_Denoise(Spectra,h,its)
% SURE_Denoise is a function that performs denoising based on Stein’s 
% Unbiased RiskEstimate. The input is the noisy signal (Spectra) 
% and the output is the denoised Spectra (denSpectra)

% Compute the size of the spectra
[A,B] = size(Spectra);

% We can use the DiscreteWavelets
% Toolbox command |ChopVector| to remove an elements from the end of
% Spectra to create vectors with a length divisible by 64. We need this to
% utilize multiple iterations of Daubechies D4 filter
for i = 1:B
    Spectra2(:,i) = ChopVector(Spectra(:,i),6);
end
[A,B] = size(Spectra2);

% Apply its iterations of the D4 filter on each column vector
Spectra3 = WT1D(Spectra2(:,2),h,its);
Spectra4 = Spectra3;

% Plot the wavelet density plot
%figure;
%WaveletVectorPlot(Spectra3,its);
%title('The wavelet transform before thresholding');

% Compute the risk value using DonohoSure in the Discrete Wavelets Toolbox
% Keep in mind this is done with only the highpass portion (lower part) of
% the transformed spectra
risk(1) = DonohoSure(Spectra3(A/2+1:A));
for i = 2:its
    risk(i) = DonohoSure(Spectra3(A/(2^i)+1:A/(2^(i-1))));
end

% At each scale, set all wavelet coefficients to 0 that are less than the
% risk value
for i = A/2+1:A
    if abs(Spectra3(i)) < risk(1)
        Spectra3(i) = 0;
    end
end
for j = 2:its
    for i = A/(2^j)+1:A/(2^(j-1))
        if abs(Spectra3(i)) < risk(j)
            Spectra3(i) = 0;
        end
    end
end

% Plot the wavelet density plot after thresholding
%figure;
%WaveletVectorPlot(Spectra3,its);
%title('The wavelet transform after thresholding');

% Perform the inverse wavelet transform
invec = IWT1D(Spectra3,h,its);

% Concatenate with the original wavelength
denSpectra = [Spectra2(:,1), invec];

end
