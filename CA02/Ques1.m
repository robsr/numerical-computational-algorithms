str=input('Enter the input filename including extension :','s');
fid = fopen(str,'r');
n=fscanf(fid,'%f',[1 1]);
A=fscanf(fid,'%f',[n+1 n]);
A=A';
fclose(fid);

choice=input('Enter 1 for Gauss Elimination(without pivoting)\nEnter 2 for Gauss Elimination(with pivoting)\nEnter 3 for Gauss Elimination(with pivoting and scaling)\nEnter 4 for LU decomposition using GE(without pivoting)\nEnter 5 for LU decomposition using GE(with pivoting)\nEnter 6 for LU decomposition using Crout Method(without pivoting)\nEnter 7 for using Cholesky decomposition(for symmetric positive definite matrix)\n');


if(choice==1)
%% Gauss Elimination(without pivoting)
    
% forward substitution
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

fprintf('Please see output_%d.txt for answers\n',choice);

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

fclose(f);
   
elseif(choice==2)
%% Gauss Elimination(with pivoting)

B=A(1:n,1:n);

% pivoting
for i = 1:n-1
    [maxc , m]= max(A(i:n,i));
    m=m+i-1;
    A([i m],:)=A([m i],:);  
end

PiA=A(1:n,1:n);   % Pivoted Matrix of coefficients only

% forward substitution
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

P = PiA*inv(B); % calculating Permutation Matrix

fprintf('Please see output_%d.txt for answers\n',choice);

% printing output

f = fopen('output_2.txt','wt');
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

fprintf(f,'\nU = \n');

for i = 1:n
    for j = 1:n
        fprintf(f,'%f ',A(i,j));
    end
    fprintf(f,'\n');
end

fclose(f);



elseif(choice==3)
%% Gauss Elimination(with pivoting and scaling)

B=A(1:n,1:n);

% Scaling
C=A;
for i=1:n
    maxr=max(C(i,:));
    C(i,:)=C(i,:)./maxr;  % making scaled matrix
end

% pivoting
for i = 1:n-1
    [maxc , m]= max(C(i:n,i)); % Pivoting according to scaled matrix
    m=m+i-1;
    A([i m],:)=A([m i],:);  
end

PiA=A(1:n,1:n);   % Pivoted Matrix of coefficients only

% forward substitution
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

P = PiA*inv(B); % calculating Permutation Matrix

fprintf('Please see output_%d.txt for answers\n',choice);

% printing output

f = fopen('output_3.txt','wt');
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

fprintf(f,'\nU = \n');

for i = 1:n
    for j = 1:n
        fprintf(f,'%f ',A(i,j));
    end
    fprintf(f,'\n');
end

fclose(f);
 
elseif(choice==4)
%% LU decomposition using GE(without pivoting)

B_real=A(:,n+1); % storing real constants B
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

fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_4.txt','wt');
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
        fprintf(f,'%f ',U(i,j));
    end
    fprintf(f,'\n');
end

fclose(f);
    
elseif(choice==5)
%% LU decomposition using GE(with pivoting)

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
    
elseif(choice==6)
%% 6 for LU decomposition using Crout Method(without pivoting)

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
    
elseif(choice==7)
%% 7 for using Cholesky decomposition(for symmetric positive definite matrix)

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
  
%%
else
    fprintf('Please enter a valid number and again run "Ques1.m"\n');
end
