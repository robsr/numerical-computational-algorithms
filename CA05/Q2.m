str = input('Enter the Required Function f(x,y):','s');
f=inline(str,'x','y');
x0=input('Enter Initial Value x0: ');
y0=input('Enter Initial Value y0: ');
xf=input('Enter Final Value xf: ');
h=input('Enter Interval Size h: ');
h_max=input('Enter Maximum Interval Size hmax: ');
alpha=input('Enter alpha for refining h: ');
tol=input('Enter Maximum Tolerance: ');


choice=input('\n\nEnter 1 for Euler Method\nEnter 2 for Mid-Point Method\nEnter 3 for 4th order RK Method\nEnter 4 for RK45 Method\n');

%% Euler's Method
if choice==1
    x(1)=x0;
    y(1)=y0;
    i=1;
    
    while x(i) <= xf
        y(i+1) = y(i)+h*f(x(i),y(i));
        x(i+1) = x(i)+h;
        i=i+1;
    end
    
    y_euler=y;
    plot(x,y_euler,'-rs');xlabel('x');ylabel('y');legend('Euler Method');grid on;hold on;
    
    % printing output
    f = fopen('out_euler.txt','wt');
   fprintf(f,'x         y_euler\n');
   
   for j = 1:length(x)
        fprintf(f,'%f  %f\n',x(j),y_euler(j));
   end
   fclose(f);
    
%% Mid Point Method
elseif choice==2
    
x(1)=x0;
y(1)=y0;
i=1;
a2=1; a1=0; p1=1/2; q11=1/2; % conditions for Mid point Method

while x(i) <= xf
        k1=f(x(i),y(i));
        k2 = f( (x(i)+h/2),(y(i)+q11*k1*h) );
        y(i+1) = y(i)+ (a1*k1 + a2*k2)*h;
        x(i+1) = x(i)+h;
        i=i+1;
end

y_midpoint=y;
plot(x,y_midpoint,'-gd');xlabel('x');ylabel('y');legend('Midpoint Method');grid on;hold on;

% printing output
    f = fopen('out_midpoint.txt','wt');
   fprintf(f,'x         y_midpoint\n');
   
   for j = 1:length(x)
        fprintf(f,'%f  %f\n',x(j),y_midpoint(j));
   end
   fclose(f);
    
%% 4th order Runge Kutta 
elseif choice==3

x(1)=x0;
y(1)=y0;
i=1;

while x(i) <= xf
        k1=f(x(i),y(i));
        k2 = f( (x(i)+h/2),(y(i)+0.5*k1*h) );
        k3 = f( (x(i)+h/2),(y(i)+0.5*k2*h) );
        k4 = f(x(i)+h, y(i)+k3*h);
        y(i+1) = y(i)+ 1/6*(k1 + 2*k2 + 2*k3 + k4)*h;
        x(i+1) = x(i)+h;
        i=i+1;
end

y_RK4=y;
plot(x,y_RK4,'-b+');xlabel('x');ylabel('y');legend('RK4 Method');grid on;hold on;

% printing output
    f = fopen('out_RK4.txt','wt');
   fprintf(f,'x         y_RK4\n');
   
   for j = 1:length(x)
        fprintf(f,'%f  %f\n',x(j),y_RK4(j));
   end
   fclose(f);

%% RK45 Method    
elseif choice==4
clear x y;

x(1)=x0;
y4(1)=y0;
y5(1)=y0;
y(1)=y0;
i=1;
h=2;

while x(i) < xf
        k1 = f(x(i),y(i));
        k2 = f(x(i)+h/5 , y(i)+0.2*k1*h);
        k3 = f(x(i)+(3/10)*h , y(i) + (3/40)*k1*h + (9/40)*k2*h);
        k4 = f(x(i)+(3/5)*h , y(i) + (3/10)*k1*h - (9/10)*k2*h + (6/5)*k3*h);
        k5 = f(x(i)+h , y(i) - (11/54)*k1*h + (5/2)*k2*h - (70/27)*k3*h + (35/27)*k4*h);
        k6 = f(x(i)+(7/8)*h , y(i) + (1631/55296)*k1*h + (175/512)*k2*h + (575/13824)*k3*h + (44275/110592)*k4*h + (253/4096)*k5*h);
        
        y4(i+1) = y(i)+ ((37/378)*k1 + (250/621)*k3 + (125/594)*k4 + (512/1771)*k6)*h;
        y5(i+1) = y(i)+ ((2825/27648)*k1 + (18575/48384)*k3 + (13525/55296)*k4 + (277/14336)*k5 + (1/4)*k6)*h;
        
        Ea = abs(y5(i+1)-y4(i+1));
        
        if Ea < tol
           
            y(i+1)=y5(i+1);
            x(i+1) = x(i)+h;
            i=i+1;
            h=2;
        else
            h = h*( (tol/Ea)^alpha );
            x(i+1) = x(i)+h;
        end
end

x(i)=4; h=4-x(i-1);

        k1 = f(x(i-1),y(i-1));
        k2 = f(x(i-1)+h/5 , y(i-1)+0.2*k1*h);
        k3 = f(x(i-1)+(3/10)*h , y(i-1) + (3/40)*k1*h + (9/40)*k2*h);
        k4 = f(x(i-1)+(3/5)*h , y(i-1) + (3/10)*k1*h - (9/10)*k2*h + (6/5)*k3*h);
        k5 = f(x(i-1)+h , y(i-1) - (11/54)*k1*h + (5/2)*k2*h - (70/27)*k3*h + (35/27)*k4*h);
        k6 = f(x(i-1)+(7/8)*h , y(i-1) + (1631/55296)*k1*h + (175/512)*k2*h + (575/13824)*k3*h + (44275/110592)*k4*h + (253/4096)*k5*h);
        
        y5(i) = y(i-1)+ ((2825/27648)*k1 + (18575/48384)*k3 + (13525/55296)*k4 + (277/14336)*k5 + (1/4)*k6)*h;
y(i)=y5(i);

y_RK45 = y;

plot(x,y_RK45,'-ko');xlabel('x');ylabel('y');legend('RK45 Method');grid on;hold on;

% printing output
    f = fopen('out_RK45.txt','wt');
   fprintf(f,'x         y_RK45\n');
   
   for j = 1:length(x)
        fprintf(f,'%f  %f\n',x(j),y_RK45(j));
   end
   fclose(f);

%%    
else
    
    printf('Enter a Valid Choice');
    
end