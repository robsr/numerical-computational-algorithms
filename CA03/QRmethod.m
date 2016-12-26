choice=4;
n=3;
A=[8 -1 -1;
    -1 4 -2;
    -1 -2 10];
toll=0.001;

A0=A;
iter=1;
A1=Rd(A0)*ON(A0);

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