function slope = secant(varargin)
% gives a normalized secant between present and last point.
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