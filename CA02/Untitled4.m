n=3;
A = [1 5 3 18 ;2 4 7 10; 4 6 2 50];
% forward substitution

B=A(1:n,1:n);

% Scaling
C=A;
for i=1:n
    maxr=max(C(i,:));
    C(i,:)=C(i,:)./maxr;  % making scaled matrix
end

% pivoting
for i = 1:n-1
    [maxc , m]= max(C(i:n,i)); %pivoting according to scale matrix
    m=m+i-1;
    A([i m],:)=A([m i],:);  
end

PiA=A(1:n,1:n);

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

P = PiA*inv(B);

% printing output

f = fopen('output_1.txt','wt');
fprintf(f,'X = \n');

for j = 1:n
        fprintf(f,'%f\n',x(j));
end

fprintf(f,'\nU = \n');

for i = 1:n
    for j = 1:n
        fprintf(f,'%f ',A(i,j));
    end
    fprintf(f,'\n');
end

fprintf(f,'\nP = \n');

for i = 1:n
    for j = 1:n
        fprintf(f,'%f ',P(i,j));
    end
    fprintf(f,'\n');
end
