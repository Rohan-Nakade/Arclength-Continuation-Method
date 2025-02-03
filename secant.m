% This function gives a normalized secant (direction) between present and last point.
% The slope calculated depends on the number of inputs:
%   If an array is given then slope is calculated using last two points.
%   If two points are given as input then the slope of these two points are
%   given.

function slope = secant(varargin)

if nargin == 2
    final = varargin{2};
    initial = varargin{1};
elseif nargin == 1
    final = varargin{1}(:,end);
    initial = varargin{1}(:,end-1);
end
grad = final-initial;
slope = grad./norm(grad);

end