function [a_grid] = make_grid(a_min,a_max,n,curv,method)

%{
INPUTS
a_min, a_max: lower and upper bound of the grid
n:            number of grid points
curv: curvature parameter.
      Method 1 If curv=1, the grid is equispaced
      If curv>1, then there are more points closer to a_min (left-skewed)
      Method 2: curv has a different interpretation!!
method: 1 = standard, 2 = Kindermann
OUTPUT
a_grid: discrete grid (as column vector)
%}

if (method==1)
    % Standard method
    a_grid = zeros(n,1);
    for i=1:n
        a_grid(i) = a_min + ((i-1)/(n-1))^curv*(a_max-a_min);
    end
    
elseif (method==2)
    % Kindermann's method
    a_grid = nan(n,1);
    % calculate factor
    h      = (a_max-a_min)/((1+growth)^n-1);
    for i=2:n+1
        a_grid(i-1) = h*((1+growth)^(i-1)-1) + a_min;
    end
    
end

end %END FUNCTION

