function z = ncc(X, Y)
X = X(:) - mean(X(:));
Y = Y(:) - mean(Y(:));
if(norm(X) == 0 || norm(Y) == 0)
    z = nan;
    return
end
z = (X' * Y) / norm(X)/norm(Y);