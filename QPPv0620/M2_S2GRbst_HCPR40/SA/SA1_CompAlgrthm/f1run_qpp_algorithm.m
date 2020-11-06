function [seed_timepoints, QPPtemplate, p] = ...
    f1run_qpp_algorithm(functional_data, number_subjects, ...
    windowlength, seed_timepoints, corr_threshold, ...
    number_iterations_thresh1, max_iterations)

% _________________________________________________________________________
%
% Algorithm developed by Waqas Majeed - modified by BY for comparison
% _________________________________________________________________________

I_msk = functional_data;
szI = size(I_msk);
nt = szI(2)/number_subjects;
QPPtemplate = zeros(szI(1), windowlength);
for k = 1:length(seed_timepoints)
    QPPtemplate = QPPtemplate + I_msk(:, seed_timepoints(k):...
        seed_timepoints(k)+windowlength-1);
end
QPPtemplate = QPPtemplate(:) - mean(QPPtemplate(:)); 
QPPtemplate = QPPtemplate / sqrt(QPPtemplate'*QPPtemplate);
correlation = zeros(szI(2)-windowlength+1, 1); 
QPPmaxima = zeros(size(correlation));
for u = 1:number_subjects
    for k = 1:nt-windowlength+1
        im_chunk = I_msk(:,(u-1)*nt+k:(u-1)*nt+k+windowlength-1); 
        im_chunk = im_chunk(:) - mean(im_chunk(:));
        correlation((u-1)*nt+k) = im_chunk'*QPPtemplate/...
            sqrt(im_chunk'*im_chunk);
    end
    sub_t = correlation( (u-1)*nt+1:(u-1)*nt+nt-windowlength+1 );
    sl = circshift(sub_t, 1);
    sr = circshift(sub_t, -1);
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
    QPPtemplate = zeros(szI(1), windowlength);
    for k = 1:length(seed_timepoints)
        QPPtemplate = QPPtemplate + I_msk(:, seed_timepoints(k):...
            seed_timepoints(k)+windowlength-1);
    end
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
        sl = circshift(sub_t, 1);
        sr = circshift(sub_t, -1);
        QPPmaxima((u-1)*nt+1:(u-1)*nt+nt-windowlength+1)=...
            ((sub_t>sl)&(sub_t>sr)).*sub_t;
        QPPmaxima((u-1)*nt+1) = 0; 
        QPPmaxima((u-1)*nt+nt-windowlength+1) = 0; 
    end
    correlation = correlation(:);
    if (ncc(t_old, correlation) > .9999) || (ncc(t_old_old, correlation)...
            > .9999)|| (ncc(t_old_old_old, correlation) > .9999)
        break;
    end
    t_old_old_old = t_old_old;
    t_old_old = t_old;
    t_old = correlation;
    p = p+1;
end
QPPtemplate = zeros(szI(1), windowlength);
for k = 1:length(seed_timepoints)
    QPPtemplate = QPPtemplate + functional_data(:, ...
        seed_timepoints(k):seed_timepoints(k)+windowlength-1);
end
QPPtemplate = QPPtemplate / length(seed_timepoints);
function z = ncc(X, Y)
X = X(:) - mean(X(:));
Y = Y(:) - mean(Y(:));
if(norm(X) == 0 || norm(Y) == 0)
    z = nan;
    return;
end
z = (X' * Y) / norm(X)/norm(Y);