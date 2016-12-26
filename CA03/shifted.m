choice=3;
n=3;
s=8;
max_it=100;
toll=0.001;
A=[8 -1 -1;
    -1 4 -2;
    -1 -2 10];

A=inv(A-s*eye(n));

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
eig=s+1/eig;


fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_3.txt','wt');
fprintf(f,'Shifted Power Method\n\n');
fprintf(f,'Eigen Value\n%f\n\n',eig);
fprintf(f,'Eigen Vector\n');

for j = 1:n
        fprintf(f,'%f\n',eig_vect(j));
end

fprintf(f,'\nIterations\n%d\n',iter);

fclose(f);
