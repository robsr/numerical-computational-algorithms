choice=6;
n=3;
A = [1 5 3 18 ;2 4 7 10; 4 6 2 50];

B_real=A(:,n+1); % storing real constants B

L=zeros(n);    % set L to be a zero matrix of n*n
U=zeros(n);    % set U to be a zero matrix of n*n

for i=1:n
    L(i,1)=A(i,1); % filling first column of L
end

for i=1:n          
    U(i,i)=1;              % diagonal elements = 1 
    U(1,i)=A(1,i)/L(1,1);  % filling first row of U
end

for j=2:n-1                % calculating remaining elements of L and U
    for i=j:n
        sum=0;
        for k=1:j-1
            sum=sum+L(i,k)*U(k,j);
        end
        L(i,j)=A(i,j)-sum;
    end
    for k=j+1:n
        sum=0;
        for i=1:j-1
            sum=sum+L(j,i)*U(i,k);
        end
        U(j,k)=( A(j,k)-sum )/L(j,j);
    end 
end

sum0=0;                 % calculating L(n,n) seperately
for k=1:n-1
    sum0=sum0+L(n,k)*U(k,n);
end
L(n,n)=A(n,n)-sum0;

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
f = fopen('output_6.txt','wt');
fprintf(f,'X = \n');

for j = 1:n
        fprintf(f,'%f\n',x(j));
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
        fprintf(f,'%f ',U_coff(i,j));
    end
    fprintf(f,'\n');
end

fclose(f);