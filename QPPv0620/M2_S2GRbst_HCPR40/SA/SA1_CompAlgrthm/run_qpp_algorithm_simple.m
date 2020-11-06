function [correlation, QPPtemplate, QPPtemplate2, QPPmaxima] = ...
    run_qpp_algorithm_simple(functional_data, mask, number_subjects, ...
    windowlength, windowlength2, seed_timepoints, corr_threshold, ...
    number_iterations_thresh1, max_iterations)

% _________________________________________________________________________
%
% This function runs a spatiotemporal pattern-finding algorithm that
% identifies quasi-periodic patterns. 
%
% Ouputs 
%
% correlation - Vector that shows sliding window correlation of QPP
% template with the functional scan inputted into algorithm. The vector is
% offset by 1/2 of the window length 
% QPPtemplate - Template of QPP 
% QPPtemplate2 - Template of QPP with longer length
% QPPmaxima - Vector that shows the maxima in the correlation vector
%
% Inputs
%
% functional_data - 2D + time functional data that is being analyzed 
% mask - Mask of data within which analysis is to be performed 
% number_subjects - The number of subjects that are concatenated in the
% entered functional data
% windowlength - Desired length of QPPtemplate
% windowlength2 - Desired length of QPPtemplate2
% seed_timepoints - Desired starting time for algorithm
% corr_threshold, number_iterations_thresh1, and max_iterations are not
% required.
%
% Algorithm developed by Waqas Majeed
% _________________________________________________________________________


if nargin < 7
    corr_threshold = [0.1 0.2];
end
if nargin < 8
    number_iterations_thresh1 = 3;
end

if nargin < 9
    max_iterations = 20;
end
szI = size(functional_data);
nt = szI(end)/number_subjects;
if round(nt)-nt > 1e-10
    error('data from multiple subjects must have the same time points ');
end
mask = (sum(abs(functional_data), 3) > 0).*mask;
ind = mask(:) > 0; %autofixed
I_msk = reshape(functional_data, szI(1)*szI(2), szI(3));
I_msk = I_msk(ind, :);
szI_msk = size(I_msk);
QPPtemplate = zeros(szI_msk(1), windowlength);
for k = 1:length(seed_timepoints)
    QPPtemplate = QPPtemplate + I_msk(:, seed_timepoints(k):...
        seed_timepoints(k)+windowlength-1);
end
QPPtemplate = QPPtemplate(:) - mean(QPPtemplate(:)); 
QPPtemplate = QPPtemplate / sqrt(QPPtemplate'*QPPtemplate);
correlation = zeros(szI_msk(2)-windowlength+1, 1); 
QPPmaxima = zeros(size(correlation));
for u = 1:number_subjects
    for k = 1:nt-windowlength+1
        im_chunk = I_msk(:,(u-1)*nt+k:(u-1)*nt+k+windowlength-1); 
        im_chunk = im_chunk(:) - mean(im_chunk(:));
        correlation((u-1)*nt+k) = im_chunk'*QPPtemplate/...
            sqrt(im_chunk'*im_chunk);
    end
    sub_t = correlation( (u-1)*nt+1:(u-1)*nt+nt-windowlength+1 );
    sl = circshift(sub_t, [1]);
    sr = circshift(sub_t, [-1]);
    QPPmaxima((u-1)*nt+1:(u-1)*nt+nt-windowlength+1) = ((sub_t>sl)&...
        (sub_t>sr)).*sub_t;
    QPPmaxima((u-1)*nt+1) = 0; QPPmaxima((u-1)*nt+nt-windowlength+1) = 0; 
end
correlation = correlation(:);
t_old = correlation;
t_old_old = correlation;
t_old_old_old = correlation;
p = 1;
while p<=max_iterations
    correlation = smooth(correlation);
    if(p <= number_iterations_thresh1)
        thresh = corr_threshold(1);
    else
        thresh = corr_threshold(2);
    end
    seed_timepoints = find((QPPmaxima) > thresh);
    
    if (length(seed_timepoints) == 1)
        break;
    end
    QPPtemplate = zeros(szI_msk(1), windowlength);
    for k = 1:length(seed_timepoints)
        QPPtemplate = QPPtemplate + I_msk(:, seed_timepoints(k):...
            seed_timepoints(k)+windowlength-1);
    end;
    QPPtemplate = QPPtemplate(:) - mean(QPPtemplate(:)); QPPtemplate = ...
        QPPtemplate / sqrt(QPPtemplate'*QPPtemplate);
    for u = 1:number_subjects
        for k = 1:nt-windowlength+1
            im_chunk = I_msk(:,(u-1)*nt+k:(u-1)*nt+k+windowlength-1); 
            im_chunk = im_chunk(:) - mean(im_chunk(:));
            correlation( (u-1)*nt+k ) = im_chunk'*QPPtemplate / ...
                sqrt(im_chunk'*im_chunk);
        end
        sub_t = correlation( (u-1)*nt+1:(u-1)*nt+nt-windowlength+1 );
        sl = circshift(sub_t, [1]);
        sr = circshift(sub_t, [-1]);
        QPPmaxima((u-1)*nt+1:(u-1)*nt+nt-windowlength+1)=...
            ((sub_t>sl)&(sub_t>sr)).*sub_t;
        QPPmaxima((u-1)*nt+1) = 0; 
        QPPmaxima((u-1)*nt+nt-windowlength+1) = 0; 
    end
    correlation = correlation(:);
    if (ncc(t_old, correlation) > .9999) || (ncc(t_old_old, correlation)...
            > .9999)|| (ncc(t_old_old_old, correlation) > .9999);
        break;
    end
    t_old_old_old = t_old_old;
    t_old_old = t_old;
    t_old = correlation;
    p = p+1;
end
QPPtemplate = zeros(szI(1), szI(2), windowlength);
for k = 1:length(seed_timepoints)
    QPPtemplate = QPPtemplate + functional_data(:, :, ...
        seed_timepoints(k):seed_timepoints(k)+windowlength-1);
end
QPPtemplate = QPPtemplate / length(seed_timepoints);
strts2 = seed_timepoints;
ind =  (round(strts2-windowlength2/2+windowlength/2+1)) < 1 ; 
strts2(ind) = [];
ind =  (round(strts2+windowlength2/2+windowlength/2) > szI(3)) ; 
strts2(ind) = [];
QPPtemplate2 = zeros(szI(1), szI(2), windowlength2);
for k = 1:length(strts2)
    QPPtemplate2 = QPPtemplate2 + functional_data(:, :, round(strts2(k)-...
        windowlength2/2+windowlength/2+1):round(strts2(k)+windowlength2/...
        2+windowlength/2));
end;
QPPtemplate2 = QPPtemplate2 / length(strts2);
function z = ncc(X, Y)
X = X(:) - mean(X(:));
Y = Y(:) - mean(Y(:));
if(norm(X) == 0 || norm(Y) == 0)
    z = nan;
    return;
end
z = (X' * Y) / norm(X)/norm(Y);