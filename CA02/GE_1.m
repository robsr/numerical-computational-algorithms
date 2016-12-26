n=3;
A=[4.0 2.0 0.0 10.0;
2.0 4.0 1.0 11.5;
0.0 1.0 5.0 4.5];
% forward substitution
for i=1:n
    if A(i,1)
    
end


for i = 1:n-1
   for j = i+1:n
       l(j,i) = A(j,i)/A(i,i);
       for k = i:n+1
           A(j,k) = A(j,k)-l(j,i)*A(i,k);
       end
   end
end
% backward elimination
x(n)=A(n,n+1)/A(n,n);

for i=1:n-1
    sum=0;
    for j=1:i
        sum=sum+A(n-i,n+1-j)*x(n+1-j);
    end
    x(n-i)=(A(n-i,n+1)-sum)/A(n-i,n-i);
end

% printing output
f = fopen('output_1.txt','wt');
A = A(:,1:n);
fprintf(f,'X = ');
dlmwrite('output_1.txt',x,' ');

dlmwrite('output_1.txt',A,' ');
fclose(f);
