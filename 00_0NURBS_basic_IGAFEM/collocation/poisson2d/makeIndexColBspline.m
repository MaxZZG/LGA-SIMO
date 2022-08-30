function [ M_arr ] = makeIndexColBspline( knot, p, deriv_order )
% pre-evaluates the Bsplines at Greville abscissae
% INPUT: knot - knot vector
%        p - polynomial degree of BSpline
%        deriv_order - number of derivatives we are evaluating
%        sample_pts - evaluation points on the reference interval [-1,1]
% OUTPUT: pt_index - connectivity array of the form 
%           pt_index(ni, loc_index#) = global_index#
%         M_arr - array of values for b_spline of the form
%           M_arr(:,:,_global_index#) = array as produced by derBasisFuns

len = length(knot)-p-1;  %number of basis functions in the u direction
point_list = zeros(1, len);
span_list = zeros(1, len);
tol = 1e-12;

for i=1:len % for each node in the x direction    
    if (i>1) && (i<len)
        point_list(i) = sum(knot(i+1:i+p))./p;   
        for j=i+1:i+p
            if (point_list(i)<=knot(j+1)) && (point_list(i)>=knot(j)) && (knot(j+1)>knot(j))
                span_list(i) = j;          
            end
        end
    end
end
span_list(end) = len;

M_arr  = evalBspline(point_list, span_list, knot, p, deriv_order);

