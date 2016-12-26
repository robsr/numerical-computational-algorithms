str=input('Enter the input filename including extension :','s');
fid = fopen(str,'r');
data=fscanf(fid,'%f',[2 5]);
data=data';
order_data=size(data);
datano=order_data(1);
xf=fscanf(fid,'%f',[1 4]);
xf=xf';
order_xf=size(xf);
xfno=order_xf(1);
s0=fscanf(fid,'%f',[1 1]);
s1=fscanf(fid,'%f',[1 1]);

fclose(fid);

z=input('Press 1 for Linear spline\nPress 2 for Quadratic spline\nPress 3 for Natural cubic spline\nPress 4 for Not-a-knot cubic spline\nPress 5 for Periodic cubic spline\nPress 6 for Clamped cubic spline\n');
figure;
hold on;

%%
if z==1
   for i=1:xfno
       j=1;
       while(1)
          if xf(i)>=data(j,1) && xf(i)<=data(j+1,1)
              break;
          end
          j=j+1;
       end
       res(i) = (xf(i)-data(j,1))/(data(j+1,1)-data(j,1))*data(j+1,2) + (xf(i)-data(j+1,1))/(data(j,1)-data(j+1,1))*data(j,2) ;
   end
    
   plot(data(:,1),data(:,2),'r*-');
   grid on;
   xlabel('x');ylabel('y');legend('Linear Spline');
   
   % printing output
   f = fopen('output_1.txt','wt');
   fprintf(f,'Interpolated values y* at give x*\n\nLinear Spline\n');
   
   for j = 1:xfno
        fprintf(f,'%f  %f\n',xf(j),res(j));
   end
   
   fclose(f);
   disp('--------------See Output_1.txt for results------------');
   
%%
elseif z==2
      A=zeros(3*datano-3,3*datano-2);
    for i=1:datano-1
            j=3*i-2;
            
            A(2*i-1,j)=data(i,1)^2;
            A(2*i,j)=data(i+1,1)^2;
            A(2*i-1,j+1)=data(i,1);
            A(2*i,j+1)=data(i+1,1);
            A(2*i,j+2)=1;
            A(2*i-1,j+2)=1;
            A(2*i-1,3*datano-2)=data(i,2);
            A(2*i,3*datano-2)=data(i+1,2);
        
    end
    for i=2*datano-1:3*datano-4
        j= 3*(i-2*datano+2)-1;
        A(i,j)=1;
        A(i,j-1)=2*data(i-2*datano+3,1);
        A(i,j+2)=-2*data(i-2*datano+3,1);
        A(i,j+3)=-1;
    end
    
    A(3*datano-3,3)=1;
    
    Areal=A(:,1:3*datano-3);
    Breal=A(:,3*datano-2);
    
    x=inv(Areal)*Breal;
    
    for i=1:xfno
       j=1;
       while(1)
          if xf(i)>=data(j,1) && xf(i)<=data(j+1,1)
              break;
          end
          j=j+1;
       end
       res(i)= (xf(i)^2) * x(3*j-2) + (xf(i))*x(3*j-1) +x(3*j) ;
    end
   
    
    for i=1:datano-1
        j=3*i-2;
        range=linspace(data(i,1),data(i+1,1));
        y=x(j)*(range.^2)+x(j+1)*range+x(j+2);
        plot(range,y,'m-');
    end
    grid on;
    xlabel('x');ylabel('y');legend('Quadratic Spline');
    
    % printing output
   f = fopen('output_2.txt','wt');
   fprintf(f,'Interpolated values y* at give x*\n\nQuadratic Spline\n');
   
   for j = 1:xfno
        fprintf(f,'%f  %f\n',xf(j),res(j));
   end
  
   disp('--------------See Output_2.txt for results------------');
 %%
elseif z==3
 
  A=zeros(datano);  
  B=zeros(datano,1);
 for i=2:datano
    h(i-1)=data(i,1)-data(i-1,1);
    g(i-1)=( data(i,2)-data(i-1,2) )/h(i-1);
 end
 
 for i=2:datano-1
    A(i,i-1)=h(i-1);
    A(i,i)=2*( h(i-1)+h(i) );
    A(i,i+1)=h(i);
    B(i)=6*( g(i)-g(i-1) );
 end

 A(1,1)=1;A(datano,datano)=1;
 
 x=inv(A)*B;
 
 for i=2:datano
     a(i-1)=x(i)/( 6*h(i-1) );
     b(i-1)=x(i-1)/( 6*h(i-1) );
     c(i-1)=( data(i,2)/h(i-1) ) - x(i)*(h(i-1))/6;
     d(i-1)=( data(i-1,2)/h(i-1) ) - x(i-1)*(h(i-1))/6;
 end
 
 for i=1:xfno
       j=1;
       while(1)
          if xf(i)>=data(j,1) && xf(i)<=data(j+1,1)
              break;
          end
          j=j+1;
       end
       res(i) = a(j)*(xf(i)-data(j,1))^3 - b(j)*(xf(i)-data(j+1,1))^3 + c(j)*(xf(i)-data(j,1)) - d(j)*(xf(i)-data(j+1,1));
 end
   
 for i=1:datano-1
        range=linspace(data(i,1),data(i+1,1));
        y=a(i)*((range-data(i,1)).^3)-b(i)*((range-data(i+1,1)).^3)+c(i)*(range-data(i,1))-d(i)*(range-data(i+1,1));
        plot(range,y,'b-');
 end
 grid on;
 xlabel('x');ylabel('y');legend('Natural Spline');
 
 % printing output
   f = fopen('output_3.txt','wt');
   fprintf(f,'Interpolated values y* at give x*\n\nNatural Spline\n');
   
   for j = 1:xfno
        fprintf(f,'%f  %f\n',xf(j),res(j));
   end
 
 disp('--------------See Output_3.txt for results------------');
   
 %%
elseif z==4
  A=zeros(datano);  
  B=zeros(datano,1);
 for i=2:datano
    h(i-1)=data(i,1)-data(i-1,1);
    g(i-1)=( data(i,2)-data(i-1,2) )/h(i-1);
 end
 
 for i=2:datano-1
    A(i,i-1)=h(i-1);
    A(i,i)=2*( h(i-1)+h(i) );
    A(i,i+1)=h(i);
    B(i)=6*( g(i)-g(i-1) );
 end

 A(1,1)=h(2); A(1,2)=-(h(1)+h(2)); A(1,3)=h(1); A(datano,datano-2)=h(datano-1); A(datano,datano-1)=-(h(datano-1)+h(datano-2)); A(datano,datano)=h(datano-2);
 
 x=inv(A)*B;
 
 for i=2:datano
     a(i-1)=x(i)/( 6*h(i-1) );
     b(i-1)=x(i-1)/( 6*h(i-1) );
     c(i-1)=( data(i,2)/h(i-1) ) - x(i)*(h(i-1))/6;
     d(i-1)=( data(i-1,2)/h(i-1) ) - x(i-1)*(h(i-1))/6;
 end
 
 for i=1:xfno
       j=1;
       while(1)
          if xf(i)>=data(j,1) && xf(i)<=data(j+1,1)
              break;
          end
          j=j+1;
       end
       res(i) = a(j)*(xf(i)-data(j,1))^3 - b(j)*(xf(i)-data(j+1,1))^3 + c(j)*(xf(i)-data(j,1)) - d(j)*(xf(i)-data(j+1,1));
 end
   
 
 for i=1:datano-1
        range=linspace(data(i,1),data(i+1,1));
        y=a(i)*((range-data(i,1)).^3)-b(i)*((range-data(i+1,1)).^3)+c(i)*(range-data(i,1))-d(i)*(range-data(i+1,1));
        plot(range,y,'c-');
 end
 grid on;
 xlabel('x');ylabel('y');legend('Not-a-knot Spline');
 
 % printing output
   f = fopen('output_4.txt','wt');
   fprintf(f,'Interpolated values y* at give x*\n\nNot-a-knot Spline\n');
   
   for j = 1:xfno
        fprintf(f,'%f  %f\n',xf(j),res(j));
   end
   
   disp('--------------See Output_4.txt for results------------');
    
 %%
elseif z==5
    
  A=zeros(datano);  
  B=zeros(datano,1);
 for i=2:datano
    h(i-1)=data(i,1)-data(i-1,1);
    g(i-1)=( data(i,2)-data(i-1,2) )/h(i-1);
 end
 
 for i=2:datano-1
    A(i,i-1)=h(i-1);
    A(i,i)=2*( h(i-1)+h(i) );
    A(i,i+1)=h(i);
    B(i)=6*( g(i)-g(i-1) );
 end

 A(1,1)=-h(1)/3; A(1,2)=-h(1)/6; A(1,datano-1)= -h(datano-1)/6; A(1,datano)=-h(datano-1)/3;
 B(1)=g(datano-1)-g(1);
 A(datano,1)=1; A(datano,datano)=-1;
 
 x=inv(A)*B;
 
 for i=2:datano
     a(i-1)=x(i)/( 6*h(i-1) );
     b(i-1)=x(i-1)/( 6*h(i-1) );
     c(i-1)=( data(i,2)/h(i-1) ) - x(i)*(h(i-1))/6;
     d(i-1)=( data(i-1,2)/h(i-1) ) - x(i-1)*(h(i-1))/6;
 end
 
 for i=1:xfno
       j=1;
       while(1)
          if xf(i)>=data(j,1) && xf(i)<=data(j+1,1)
              break;
          end
          j=j+1;
       end
       res(i) = a(j)*(xf(i)-data(j,1))^3 - b(j)*(xf(i)-data(j+1,1))^3 + c(j)*(xf(i)-data(j,1)) - d(j)*(xf(i)-data(j+1,1));
 end
 
 
 for i=1:datano-1
        range=linspace(data(i,1),data(i+1,1));
        y=a(i)*((range-data(i,1)).^3)-b(i)*((range-data(i+1,1)).^3)+c(i)*(range-data(i,1))-d(i)*(range-data(i+1,1));
        plot(range,y,'y-');
 end
 grid on;
 xlabel('x');ylabel('y');legend('Periodic Spline');
 
 % printing output
   f = fopen('output_5.txt','wt');
   fprintf(f,'Interpolated values y* at give x*\n\nPeriodic Spline\n');
   
   for j = 1:xfno
        fprintf(f,'%f  %f\n',xf(j),res(j));
   end
   
   disp('--------------See Output_5.txt for results------------');
   
 %%
elseif z==6
    
  A=zeros(datano);  
  B=zeros(datano,1);
 for i=2:datano
    h(i-1)=data(i,1)-data(i-1,1);
    g(i-1)=( data(i,2)-data(i-1,2) )/h(i-1);
 end
 
 for i=2:datano-1
    A(i,i-1)=h(i-1);
    A(i,i)=2*( h(i-1)+h(i) );
    A(i,i+1)=h(i);
    B(i)=6*( g(i)-g(i-1) );
 end

 A(1,1)=-h(1)/3; A(1,2)=-h(1)/6;
 B(1)=s0-g(1);
 A(datano,datano)=h(datano-1)/3; A(datano,datano-1)=h(datano-1)/6;
 B(datano)=s1-g(datano-1);
 
 x=inv(A)*B;
 
 for i=2:datano
     a(i-1)=x(i)/( 6*h(i-1) );
     b(i-1)=x(i-1)/( 6*h(i-1) );
     c(i-1)=( data(i,2)/h(i-1) ) - x(i)*(h(i-1))/6;
     d(i-1)=( data(i-1,2)/h(i-1) ) - x(i-1)*(h(i-1))/6;
 end
 
 for i=1:xfno
       j=1;
       while(1)
          if xf(i)>=data(j,1) && xf(i)<=data(j+1,1)
              break;
          end
          j=j+1;
       end
       res(i) = a(j)*(xf(i)-data(j,1))^3 - b(j)*(xf(i)-data(j+1,1))^3 + c(j)*(xf(i)-data(j,1)) - d(j)*(xf(i)-data(j+1,1));
 end
   
 
 for i=1:datano-1
        range=linspace(data(i,1),data(i+1,1));
        y=a(i)*((range-data(i,1)).^3)-b(i)*((range-data(i+1,1)).^3)+c(i)*(range-data(i,1))-d(i)*(range-data(i+1,1));
        plot(range,y,'g-');
 end
 xlabel('x');ylabel('y');legend('Clamped Spline');
 grid on;
 
 % printing output
   f = fopen('output_6.txt','wt');
   fprintf(f,'Interpolated values y* at give x*\n\nClamped Spline\n');
   
   for j = 1:xfno
        fprintf(f,'%f  %f\n',xf(j),res(j));
   end
  
   disp('--------------See Output_6.txt for results------------');
end