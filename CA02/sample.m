

B=A(1:n ,1:n);
n=3;
for i = 1:n-1
    % pivoting
    [maxc , m]= max(A(i:n,i));
    A([1 m],:)=A([m 1],:);
    
   for j = i+1:n
       l(j,i) = A(j,i)/A(i,i);
       for k = i:n+1
           A(j,k) = A(j,k)-l(j,i)*A(i,k);
       end
   end
end
P = A(1:n,1:n)*inv(B);
disp(P);