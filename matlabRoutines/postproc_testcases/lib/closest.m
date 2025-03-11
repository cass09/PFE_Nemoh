function [c,x_out,err]=closest(x,x_in,n2)
n1 =length(x_in);
if nargin == 2,n2=1;end
x_out = zeros(n1,n2);
err = zeros(n1,n2);
c = zeros(n1,n2);
n = 1:length(x);
for i=1:n1
    a=n;
    for j=1:n2
        if j==1
            [err(i,j),b] = min(abs(x-x_in(i)));
            c(i,j) = b;
        else
            a = setdiff(a,c(i,j-1));
            [err(i,j),b] = min(abs(x(a)-x_in(i)));
            c(i,j)= a(b);
        end
        x_out(i,j) = x(c(i,j));
    end
end