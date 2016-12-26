str=input('Enter the input filename including extension :','s');
fid = fopen(str,'r');
n=fscanf(fid,'%f',[1 1]);
A=fscanf(fid,'%f',[n n]);
A=A';
max_it=fscanf(fid,'%d',[1 1]);
toll=fscanf(fid,'%f',[1 1]);
s=fscanf(fid,'%f',[1 1]);
fclose(fid);

choice=input('Enter 1 for Direct Power Method\nEnter 2 for Inverse Power Method\nEnter 3 for Shifted-Power Method\nEnter 4 for QR Method\n');

if choice==1
%% Direct Power Method

X=ones(1,n);
X=X';
X_new=A*X;
[eig,index] = max(abs(X_new));    % finding maximum abs value of eigen value
X = X_new./X_new(index);          % dividing X_new by maximum absolute eigen value
error=100;
pre_eig=eig;
iter=1;

while ( error>=toll && iter<=max_it)
      iter=iter+1;
      X_new=A*X;
      [eig,index] = max(abs(X_new));
      eig_vect=X;
      error=abs(eig-pre_eig)*100/eig;
      pre_eig=eig;
      X = X_new./X_new(index);
end


fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_1.txt','wt');
fprintf(f,'Direct Power Method\n\n');
fprintf(f,'Eigen Value\n%f\n\n',eig);
fprintf(f,'Eigen Vector\n');

for j = 1:n
        fprintf(f,'%f\n',eig_vect(j)/norm(eig_vect));
end

fprintf(f,'\nIterations\n%d\n',iter);

fclose(f);
    
elseif choice==2
%% Inverse Power Method

A=inv(A);                                    % taking A to be inverse of A
X=ones(1,n);
X=X';
X_new=A*X;
[eig,index] = max(abs(X_new));
X = X_new./X_new(index);
error=100;
pre_eig=eig;
iter=1;

while ( error>=toll && iter<=max_it)
      iter=iter+1;
      X_new=A*X;
      [eig,index] = max(abs(X_new));
      eig_vect=X;
      error=abs(eig-pre_eig)*100/eig;
      pre_eig=eig;
      X = X_new./X_new(index);
end
eig=1/eig;                              % determining original eigen value by reciprocating

fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_2.txt','wt');
fprintf(f,'Inverse Power Method\n\n');
fprintf(f,'Eigen Value\n%f\n\n',eig);
fprintf(f,'Eigen Vector\n');

for j = 1:n
        fprintf(f,'%f\n',eig_vect(j)/norm(eig_vect));
end

fprintf(f,'\nIterations\n%d\n',iter);

fclose(f);
    
elseif choice==3
%% Shifted-Power Method

A=[8 -1 -1;
    -1 4 -2;
    -1 -2 8];

A=inv(A-s*eye(n));                                 % shifting matrix A by s

X=ones(1,n);
X=X';
X_new=A*X;
[eig,index] = max(abs(X_new));
X = X_new./X_new(index);
error=100;
pre_eig=eig;
iter=1;

while ( error>=toll && iter<=max_it)
      iter=iter+1;
      X_new=A*X;
      [eig,index] = max(abs(X_new));
      eig_vect=X;
      error=abs(eig-pre_eig)*100/eig;
      pre_eig=eig;
      X = X_new./X_new(index);
end
eig=s+1/eig;                                    % determining original eigen value by adding s and reciprocating


fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_3.txt','wt');
fprintf(f,'Shifted Power Method\n\n');
fprintf(f,'Eigen Value\n%f\n\n',eig);
fprintf(f,'Eigen Vector\n');

for j = 1:n
        fprintf(f,'%f\n',eig_vect(j)/norm(eig_vect));
end

fprintf(f,'\nIterations\n%d\n',iter);

fclose(f);

elseif choice==4
%% QR Method

A0=A;
iter=1;
A1=Rd(A0)*ON(A0);        % calling functions Rd and ON where ON and Rd gives Q and R such that A=Q*R 
iter=0;

while 1
      iter=iter+1;
      for i=1:n
          errors(i)=abs( (A1(i,i)-A0(i,i))*100 )/A1(i,i);
          
      end
    err_max=max(errors);
    
    if err_max<toll
        break;
    end
    
    A0=A1;
    A1=Rd(A1)*ON(A1);
end

fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_4.txt','wt');
fprintf(f,'QR Method\n\n');
fprintf(f,'Eigen Value\n');
for j = 1:n
        fprintf(f,'%f\n',A1(j,j));
end
fprintf(f,'\nIterations\n%d\n',iter);

fclose(f);
  
else
%%
    fprintf('Please Enter a Valid Number\n');
    
end