function [x] = Smooth(Spectra,piece,n)
% This is a matlab function that perfroms signal smoothing using an
% averageing filter procedure.

% The input arguments are Spectra, the desired signal to be smoothed, and
% Piece, the number of vector components per filter iteration

% Compute the size of Spectra
[A,B] = size(Spectra);

x = Spectra;
Spectra = [Spectra(:,B-1) Spectra(:,B)];
% Perform signal smoothing
for i = 1:n
    smooth = ones(1,piece)/piece;
    smoothsignal = filter(smooth,1,Spectra(:,2));

    % Adjust for filter delay
    fdelay = (Spectra(piece,1)-Spectra(1,1))/2;
    Spectra = [Spectra(:,1)-fdelay, smoothsignal];
    x = [x, Spectra];
end
end

