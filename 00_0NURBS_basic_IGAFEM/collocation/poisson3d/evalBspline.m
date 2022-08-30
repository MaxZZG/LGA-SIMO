function [ M_arr ] = evalBspline(pts, span_arr, knotu, p, deriv_order)
%pre-evaluate the B-splines in 1d for a given knot vector at all the
%evaluation points
%INPUT: pts - vector or points in parameter space to be evaluated
%       span_arr - the knot span index corresponding to each entry in pts
%       knotu - knot vector
%       p - polynomial degree
%       deriv order - number of derivatives (0: only eval shape functions,
%       1 eval shape functions and 1st derivative, etc.
%OUTPUT: M_arr(i,j,k) - M_arr(:,:,k) - array of shape funcions and or
%        derivatives as returned by derBasisFuns for the kth sample point

num_pts = length(pts);
M_arr = zeros(deriv_order+1, p+1, num_pts);
for i=1:length(pts)
    M_arr(:,:,i)=dersbasisfuns3(span_arr(i), p, deriv_order, pts(i), knotu);
    %M_arr(:,:,i)=derBasisFuns(span_arr(i), p, deriv_order, pts(i), knotu);
end

