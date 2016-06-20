function y = rd_wmean(x,w,dim)
%
% function y = rd_wmean(x,w,dim)
%
% Take a weighted mean of x.
%
% x is a matrix (up to 4D)
% w is a vector of weights
% dim is the dimension of x over which to take the weighted mean. must be 1
% or 2. default is 1.

if nargin<3
    dim = 1;
end
if dim~=1 && dim~=2
    error('dim must be 1 or 2')
end
if ~isvector(w)
    error('w must be a vector (1xn)')
end
if size(x,dim)~=length(w)
    error('the length of w must be equal to the size of x in dimension dim. default dim = 1')
end

% make w a scaled column vector
w = w/sum(w);
if size(w,1)==1
    w = w';
end

d = size(x);
ndims = numel(d);

switch ndims
    case {1,2}
        if dim==1
            y(1,:) = x'*w;
        else
            y(:,1) = x*w;
        end
    case 3
        for i = 1:d(3)
            x0 = x(:,:,i);
            if dim==1
                y(1,:,i) = x0'*w;
            else
                y(:,1,i) = x0*w;
            end
        end
    case 4
        for i = 1:d(3)
            for j = 1:d(4)
                x0 = x(:,:,i,j);
                if dim==1
                    y(1,:,i,j) = x0'*w;
                else
                    y(:,1,i,j) = x0*w;
                end
            end
        end
    otherwise
        error('x cannot have more than 4 dimensions')
end
            