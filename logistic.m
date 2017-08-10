function out = logistic(varargin)

% get params
varargin{4} =[];
if ~isempty(varargin{1}), x0=varargin{1};else x0=rand;end
if ~isempty(varargin{2}), a=varargin{2};else a=3.5;end
if ~isempty(varargin{3}), n=varargin{3};else n=100;end

out = zeros( n, 1 );
out(1) = x0;
for i = 1:n-1
	out(i+1) = a*out(i)*(1-out(i));
end