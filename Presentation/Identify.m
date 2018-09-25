function [ ] = Identify(peaks)
% This is a Matlab frunction that identifies the functional groups present
% given a list of peak values.

% peaks should be a matrix with two columns, the first column corresponding
% to the wavelenth, and the second column corresponding to the percent
% transmittence

% Compute the size of peaks
[A,B] = size(peaks);

% For each peak, find a functional group
for i = 1:A
    if peaks(i,2) < 95
        
    % Carboxylic Acids
    if peaks(i,1) <= 3200 && peaks(i,1) >= 2500
        disp(sprintf('Carboxylic Acid'));
    end
    
    % Alcohols
    if peaks(i,1) <= 3650 && peaks(i,1) >= 3200
        disp(sprintf('Alcohol'));
    end

    % Carbon-Hydrogen Bonds
    if peaks(i,1) <= 3300 && peaks(i,1) >= 2700
        disp(sprintf('Methyl Group'));
    end
    
    % Carbon-Nitrogen Bonds and Ether
    if peaks(i,1) <= 1250 && peaks(i,1) >= 1020
        disp(sprintf('Ether or C-N Bond'));
    end
    
    % Carobonyl
    if peaks(i,1) <= 1780 && peaks(i,1) >= 1640
        disp(sprintf('Carbonyl'));
    end
    
    % Aromatic Ring
    if peaks(i,1) <= 1500 && peaks(i,1) >= 1430
        disp(sprintf('Aromatic Ring'));
    end
    
    % Carbon-Nitrogen double bond
    if peaks(i,1) <= 1650 && peaks(i,1) >= 1550
        disp(sprintf('Carbon-Nitrogen Double Bond'));
    end
 
    % Alkene
    if peaks(i,1) <= 1680 && peaks(i,1) >= 1600
        disp(sprintf('Alkene'));
    end
    
    % Alkyne
    if peaks(i,1) <= 2260 && peaks(i,1) >= 2100
        disp(sprintf('Alkyne'));
    end
    
    % Carbon-hydrogen bond
    if peaks(i,1) <= 770 && peaks(i,1) >= 630
        disp(sprintf('Carbon-Hydrogen Bond'));
    end
    end
end
end

