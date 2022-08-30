function [ pt_index, M_arr ] = makeIndexBspline( knot, p, deriv_order, sample_pts )
% pre-evaluates the Bsplines at given sample_pts (on the reference interval
% [-1, 1]) on each knot span
% INPUT: knot - knot vector
%        p - polynomial degree of BSpline
%        deriv_order - number of derivatives we are evaluating
%        sample_pts - evaluation points on the reference interval [-1,1]
% OUTPUT: pt_index - connectivity array of the form 
%           pt_index(ni, loc_index#) = global_index#
%         M_arr - array of values for b_spline of the form
%           M_arr(:,:,_global_index#) = array as produced by derBasisFuns

tol = 1e-12;
span_counter = 0;
num_points = length(sample_pts);
num_uniq_span = length(unique(knot))-1;
num_span = length(knot)-1;
list_size = num_uniq_span*num_points;
point_list = zeros(1, list_size);
span_list = zeros(1, list_size);
pt_index = zeros(num_span,num_points);

point_counter = 0;
for i=1:num_span
    if abs(knot(i+1)-knot(i))>tol
        span_counter = span_counter + 1;
        %map the points on the reference interval to the parameter space
        for j=1:num_points
            point_counter = point_counter + 1;
            u = ((knot(i+1)-knot(i))*sample_pts(j) +knot(i+1) + knot(i))/2;
            point_list(point_counter) = u;
            span_list(point_counter) = i;
            pt_index(i,j) = point_counter;
            
        end
    end

end
M_arr  = evalBspline(point_list, span_list, knot, p, deriv_order);

