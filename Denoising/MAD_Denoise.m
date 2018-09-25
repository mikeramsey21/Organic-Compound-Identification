function [denSpectra] = MAD_Denoise(Spectra,h,its)
% MAD_Denoise is a function that performs denoising based on the Median
% absolute deviation method. The input is the noisy signal (Spectra) 
% and the output is the denoised Spectra (denSpectra)

% h is the filter type and its is the number of iterations

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

% Apply its interations of the D4 filter on each column vector
Spectra3 = WT1D(Spectra2(:,2),h,its);
Spectra4 = Spectra3;

% Calculate the median absolute deviation for each scale
s(1) = median(abs(Spectra3(A/2+1:A)));
for i = 2:its
    s(i) = median(abs(Spectra3(A/(2^i)+1:A/(2^(i-1)))))/.6745;
end

% Divide out by the scale
Spectra3(A/2+1:A) = Spectra3(A/2+1:A) ./ s(1);
for i = 2:its
    Spectra3(A/(2^i)+1:A/(2^(i-1))) = Spectra3(A/(2^i)+1:A/(2^(i-1))) ./ s(i);
end

% Calculate the thresholds at each scale
t(1) = (2 * log(length(Spectra3(A/2+1:A))))^.5;
for i = 2:its
    t(i) = (2 * log(length(Spectra3(A/(2^i)+1:A/(2^(i-1))))))^.5;
end

% Perform soft thresholding
for i = A/2+1:A
    if abs(Spectra3(i)) > t(1)
        if Spectra3(i) < 0
            Spectra4(i) = -(abs(Spectra4(i)) - t(1));
        elseif Spectra3(i) > 0 
            Spectra4(i) = abs(Spectra4(i)) - t(1);
        end
    else
        Spectra4(i) = 0;
    end
end
for j = 2:its
    for i = 1:length(Spectra3(A/(2^j)+1:A/(2^(j-1))))
        if abs(Spectra3(i)) > t(j)
            if Spectra3(i) < 0
                Spectra4(i) = -(abs(Spectra4(i)) - t(j));
            elseif Spectra3(i) > 0 
                Spectra4(i) = abs(Spectra4(i)) - t(j);
            end
        else
            Spectra4(i) = 0;
        end
    end
end

% Perform the invers wavelet transform
invec = IWT1D(Spectra4,h,its);

% Concatenate with the original wavelength
denSpectra = [Spectra2(:,1), invec];
end

