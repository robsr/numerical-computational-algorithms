function Q1=Q1()

str = input('Enter the Required Function:','s');
global f;
f=inline(str);
a=input('Enter lower limit a: ');
b=input('Enter upper limit b: ');
tol=input('Enter maximum allowable approximate error for each sub interval: ');

global count;
count=3;

c=(a+b)/2;
fa=f(a);
fc=f(c);
fb=f(b);
global arr;
arr=[a b c];

answ=integral(a,b,c,fa,fb,fc,tol);

fprintf('\nI = %f\n',answ);
fprintf('n = %f\n',count-1);
x=sort(arr);

for i=1:count
    y(i)=f(x(i));
end

plot(x,y,'*b-');xlabel('x');ylabel('f(x)');
grid on;

end
function integ= integral(a,b,c,fa,fb,fc,tol)
    d=(a+c)/2;
    e=(b+c)/2;
    global f;
    global arr;
    global count;
    count=count+2;
    fd=f(d);
    fe=f(e);
    arr=[arr d e];
    
    I1=(b-a)*(fa+4*fc+fb)/6;
    I2=(b-a)*(fa+4*fd+2*fc+4*fe+fb)/12;
    
    if abs(I2-I1)<=tol
        integ=I2+(I2-I1)/15;
    else
        Ia=integral(a,c,d,fa,fc,fd,tol);
        Ib=integral(c,b,e,fc,fb,fe,tol);
        integ=Ia+Ib;
    end
   end