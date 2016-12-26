choice=5;
n=3;
a = [1 5 3 18 ;2 4 7 10; 4 6 2 50];

for k=1:n
 l(k,1)=a(k,1);
 end
 
 for j=1:n
 u(1,j)=a(1,j)/l(1,1);
 u(j,j)=1;
 end
 
 for r=2:n
 
    for t=r:n
 
    sum=0;
 
        for k=1:r-1
 
            sum=l(t,k)*u(k,r);
 
        end
    l(t,r)=a(t,r)-sum;
 
    end
 
    for t=r+1:n
 
        sum=0;
    
        for k=1:r-1
 
            sum=u(k,t)*l(r,k);
 
        end
 
        u(r,t)=(a(r,t)-sum)/l(r,r);
 
    end
 
 end 
 
 disp(l);
 disp(u);
 %disp(l*u);
 
 
 k=1;
while(k<=n)

sum=0;

for j=1:k-1

sum=sum+l(k,j)*y(j);

end
y(k)=(a(k,n+1)-sum)/l(k,k);

k=k+1;
end


k=n;
while(k>0)

sum=0;

for j=k+1:n

sum = sum + u(k,j)*x(j);

end
x(k)=(y(k)-sum)/u(k,k);

k=k-1;
end
disp(y);
disp(x);