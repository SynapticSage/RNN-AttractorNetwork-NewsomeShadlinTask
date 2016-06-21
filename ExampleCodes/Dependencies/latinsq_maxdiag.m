function L = latinsq_maxdiag(N)
% This function calls the latin square function and manipulates it so that
% the maximum entries live on the diagonal.

% Obtain a latin square of size N
L = latsq(N);

% Move the maxima back to their natural diagonal
L = L(end:-1:1,1:end);

end