%%%%%%Generate Grid with n*n points(Base Stations);
%%%%%%Output the X axis and Y axis of all points.

function B = GridBaseLayout(r,n)
N = n^2;
X = zeros(n,1);
Y = zeros(n,1);
A = zeros(n,n,2);
for i = 1 : n
    X(i) = r/n*(i-1)+r/2/n;
    Y(i) = r/n*(i-1)+r/2/n;
end
for i = 1 : n
    for j = 1 : n
        A(i,j,:) = [X(i)-r/2 Y(j)-r/2];
    end
end
B = reshape(A,n^2,2);
        
