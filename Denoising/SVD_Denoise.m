function [denSpectra] = SVD_Denoise(Spectra,k)
% SVD_Denoise is a function that performs denoising based on the singular
% value Decomposition Algorithm. The input is the noisy signal (Spectra) 
% and the output is the denoised Spectra (denSpectra)

% Compute the size of the spectra
[A,B] = size(Spectra);

% Declare some modification variables
% Options for piece: 27, 81, 89
piece = 89; num = A/piece;

% We begin to perform the singular value decomposition
meth = Spectra(:,2); wave = Spectra(:,1);

% We need to create the Hankel matrix for the vector above
% We divde the signal vector into 34 chunks of length 212
for i = 1:piece
    Temp(:,i) = meth(num*(i-1)+1:num*i);
    Temp2(:,i) = wave(num*(i-1)+1:num*i);
end

% Calculate the midpoint of each piece
mid = int32(length(Temp(:,1))/2+.5);

% Compute the Hankel Matrix for each piece of signal
for i = 1:piece
    Hankel(:,:,i) = hankel(Temp(1:mid,i),Temp(mid:num,i));
end

% Compute the SVD for each Hankel Matrix
for i = 1:piece
    [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(Hankel(:,:,i));
end

% Suppose we keep the first k singular values
for i = 1:piece
    Hankel2(:,:,i) = U(:,1:k,i)*S(1:k,1:k,i)*V(:,1:k,i)';
end

% We convert the Hankel Matrices back into vectors
for i = 1:piece
    Vec(:,i) = [Hankel2(:,1,i); Hankel2(2:mid,mid,i)];
end

% Concatenate the Vec and Temp2 pieces into one vector
revspec = []; signal = [];
for i = 1:piece
    revspec = [revspec; Vec(:,i)];
    signal = [signal; Temp2(:,i)]; 
end

% Finally se concatenate the wavelength and spectra pieces
denSpectra = [signal, revspec];

end

