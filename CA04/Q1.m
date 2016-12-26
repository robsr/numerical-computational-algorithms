data=load('input1.txt');
n=input('Enter Degree of Polynomial:');
order=size(data);
N=order(1);
A=zeros(n+1,n+2);

for i=1:n+1
 for j=1:n+1;
     sum=0;
     for k=1:N
         sum=sum+data(k,1)^(i+j-2);
     end
     A(i,j)=sum;
 end
end

for i=1:n+1;
    sum=0;
    for j=1:N
        sum=sum+data(j,2)*data(j,1)^(i-1);
    end
    A(i,n+2)= sum;
end


%% Gauss Elimination(with pivoting and scaling)

% Scaling
C=A;
for i=1:n+1
    maxr=max(C(i,:));
    C(i,:)=C(i,:)./maxr;  % making scaled matrix
end

% pivoting
for i = 1:n
    [maxc , m]= max(C(i:n+1,i)); % Pivoting according to scaled matrix
    m=m+i-1;
    A([i m],:)=A([m i],:);  
end


% forward substitution
for i = 1:n
   for j = i+1:n+1
       l(j,i) = A(j,i)/A(i,i);
       for k = i:n+2
           A(j,k) = A(j,k)-l(j,i)*A(i,k);
       end
   end
end

% backward elimination
x(n+1)=A(n+1,n+2)/A(n+1,n+1);

for i=1:n
    sum=0;
    for j=1:i
        sum=sum+A(n+1-i,n+2-j)*x(n+2-j);
    end
    x(n+1-i)=(A(n+1-i,n+2)-sum)/A(n+1-i,n+1-i);
end

Yavg=mean(data(:,2));
arr0=( (data(:,2)-Yavg).^2 );

s0=0;
for i=1:N
    s0=s0+arr0(i);
end 



for i=1:N
    sum=0;
    for j=1:n+1
        sum=sum+x(j)*data(i,1)^(j-1);
    end
    y(i,1)= sum;
end


arr1= ((data(:,2)-y(:,1)).^2);

s1=0;
for i=1:N
    s1=s1+arr1(i);
end

Rsqr=1-(s1/s0);



% plotting
range=linspace(0,1);
y=zeros(size(range));

for i=1:n+1
    y=y+x(i)*range.^(i-1);
end

plot(range,y,'b-');
hold on;
plot(data(:,1),data(:,2),'*r');xlabel('x');ylabel('y');

% printing output
   f = fopen('output.txt','wt');
   if n==1
       fprintf(f,'Linear : Coefficients   ');
       for j = 1:n+1
        fprintf(f,'%f   ',x(j));
       end
       fprintf(f,'\n:R-sq = %f',Rsqr);
   
   elseif n==2
       fprintf(f,'Quadratic : Coefficients   ');
       for j = 1:n+1
        fprintf(f,'%f   ',x(j));
       end
       fprintf(f,'\n:R-sq = %f',Rsqr);
       
   elseif n==3
       fprintf(f,'Cubic : Coefficients   ');
       for j = 1:n+1
        fprintf(f,'%f   ',x(j));
       end
       fprintf(f,'\n:R-sq = %f',Rsqr);
       
   elseif n==4
       fprintf(f,'Quartic : Coefficients   ');
       for j = 1:n+1
        fprintf(f,'%f   ',x(j));
       end
       fprintf(f,'\n:R-sq = %f',Rsqr);
       
   end
      
   fclose(f);
   disp('--------------See Output.txt for results------------');













    

