choice=5;
n=3;
A = [1 5 3 18 ;2 4 7 10; 4 6 2 50];



B=A(1:n,1:n);

% pivoting
for i = 1:n-1
    [maxc , m]= max(A(i:n,i));
    m=m+i-1;
    A([i m],:)=A([m i],:);  
end
B_real=A(:,n+1); % storing real constants B
PiA=A(1:n,1:n);   % Pivoted Matrix of coefficients only

L=zeros(n);    % set L to be a zero matrix of n*n
for i=1:n
    L(i,i)=1;  % set diagonal elements to be 1
end

% forward substitution
for i = 1:n-1
   for j = i+1:n
       L(j,i) = A(j,i)/A(i,i); % filling L with rest of the elements
       for k = i:n+1
           A(j,k) = A(j,k)-L(j,i)*A(i,k);
       end
   end
end

L_coff=L;
U = A(:,1:n);  % defining U upper triangular

% calculating y from Ly=b by forward elimination

L=[L B_real]; % making augmented matrix for L

y(1)=L(1,n+1)/L(1,1);

for i=1:n-1
    sum=0;
    for j=1:i
        sum=sum+L(i+1,j)*y(j);
    end
    y(i+1)=(L(i+1,n+1)-sum)/L(i+1,i+1);
end


% calculating x from Ux=y by backward elimination

A(:,n+1)=y';  %making augmented matrix for U upper triangular

x(n)=A(n,n+1)/A(n,n);

for i=1:n-1
    sum=0;
    for j=1:i
        sum=sum+A(n-i,n+1-j)*x(n+1-j);
    end
    x(n-i)=(A(n-i,n+1)-sum)/A(n-i,n-i);
end

P = PiA*inv(B); % calculating Permutation Matrix

fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_5.txt','wt');
fprintf(f,'X = \n');

for j = 1:n
        fprintf(f,'%f\n',x(j));
end

fprintf(f,'\nPermutation Matrix = \n');

for i = 1:n
    for j = 1:n
        fprintf(f,'%f ',P(i,j));
    end
    fprintf(f,'\n');
end

fprintf(f,'\nL = \n');

for i = 1:n
    for j = 1:n
        fprintf(f,'%f ',L_coff(i,j));
    end
    fprintf(f,'\n');
end

fprintf(f,'\nU = \n');

for i = 1:n
    for j = 1:n
        fprintf(f,'%f ',U(i,j));
    end
    fprintf(f,'\n');
end

fclose(f);