choice=7;
n=3;
A=[1 5 3 18;
2 4 7 10;
4 6 2 50];

B_real=A(:,n+1); % storing real constants B

L=zeros(n);    % set L to be a zero matrix of n*n
L(1,1)=sqrt( A(1,1) );

for i=1:n
    L(i,1)=A(i,1)/L(1,1); % filling first column of L
end

for j=2:n              % calculating remaining elements of L
    sum=0;
    for k=1:j-1
        sum=sum+L(j,k)*L(j,k);
    end
    L(j,j) = sqrt( A(j,j)-sum );
    for k=j+1:n
        sum=0;
        for s=1:j-1
            sum=sum+L(j,s)*L(k,s);
        end
        L(k,j)=( A(k,j)-sum )/L(j,j);
    end 
end
U = L';               % defining Upper triangular matrix 

L_coff=L;             % storing L and U in n*n matrices
U_coff=U;

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

U=[U y'];  %making augmented matrix for U upper triangular

x(n)=U(n,n+1)/U(n,n);

for i=1:n-1
    sum=0;
    for j=1:i
        sum=sum+U(n-i,n+1-j)*x(n+1-j);
    end
    x(n-i)=(U(n-i,n+1)-sum)/U(n-i,n-i);
end

fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_7.txt','wt');
fprintf(f,'X = \n');

for j = 1:n
        fprintf(f,'%f\n',x(j));
end

fprintf(f,'\nCholesky Factor = \n');

for i = 1:n
    for j = 1:n
        fprintf(f,'%f ',L_coff(i,j));
    end
    fprintf(f,'\n');
end

fclose(f);