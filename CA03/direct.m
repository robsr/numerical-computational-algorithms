choice=1;
n=3;
max_it=100;
toll=0.001;
A=[8 -1 -1;
    -1 4 -2;
    -1 -2 10];
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


fprintf('Please see output_%d.txt for answers\n',choice);

% printing output
f = fopen('output_1.txt','wt');
fprintf(f,'Direct Power Method\n\n');
fprintf(f,'Eigen Value\n%f\n\n',eig);
fprintf(f,'Eigen Vector\n');

for j = 1:n
        fprintf(f,'%f\n',eig_vect(j));
end

fprintf(f,'\nIterations\n%d\n',iter);

fclose(f);
