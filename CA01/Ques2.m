
k=input('Press 1 for Muller Method\nPress 2 for Bairstow Method\n');

if (k==1)
%% Muller Method

str=input('Enter the given Polynomial :','s');
f=inline(str);

x0=input('Enter the first Guess Value :');x1=input('Enter the second Guess Value :');x2=input('Enter the third Guess Value :');
   maxit=input('Enter Maximum no.of Iterations :');errmax=input('Enter the Maximum Relative approximate error in(%) :');
   
   err=[];
   flag=0;
   
   
   if (f(x0)==0)
           fprintf('Root of given Equation is :%f\n',x0);
           flag=1;
   elseif(f(x1)==0)
           fprintf('Root of given Equation is :%f\n',x1);
            flag=1;
   elseif(f(x2)==0)
           fprintf('Root of given Equation is :%f\n',x2);
           flag=1;
   else
       h0=x1-x0;
       h1=x2-x1;
       f0=(f(x1)-f(x0))/(x1-x0);
       f1=(f(x2)-f(x1))/(x2-x1);
       a=(f1-f0)/(h1+h0);
       b=a*h1 + f1;
       c=f(x2);
       
       if(b>0)
           x3 = x2 + (-2*c)/(b+sqrt(b^2-4*a*c));
       else
           x3 = x2 + (-2*c)/(b-sqrt(b^2-4*a*c));
       end
       
       prex0=x0;
       prex1=x1;
       prex2=x2;
       prex3=x3;
       
        for i=1:maxit
            if (f(prex3)==0)
                fprintf('Root of given Equation is :%f',x3);
                flag=1;
                break;
            else
                h0=prex2-prex1;
                h1=prex3-prex2;
                f0=(f(prex2)-f(prex1))/(prex2-prex1);
                f1=(f(prex3)-f(prex2))/(prex3-prex2);
                a=(f1-f0)/(h1+h0);
                b=a*h1 + f1;
                c=f(prex3);
                
                if(b>0)
                       currx3 = prex3 + (-2*c)/(b+sqrt(b^2-4*a*c));
                else
                       currx3 = prex3 + (-2*c)/(b-sqrt(b^2-4*a*c));
                end
                
                err(i)=abs((currx3-prex3)/currx3)*100;
                
                prex1=prex2;
                prex2=prex3;
                prex3=currx3;
            end
          if(err(i)<errmax)
              fprintf('\nRoot of given Equation is :%f\n',currx3);
              flag=1;
              break;
          end
        end
      if flag==0
       fprintf('Solution not found\nNo. of maximum iterations provided are maybe less\nTry with more no. of iterations\n');
      end
      
   x=-10:10;
   figure;
   plot(x,f(x),'-b');grid on;xlabel('x');ylabel('f(x)');title('f(x) vs x');
   legend('f(x)');    
   ax=gca;
   ax.XAxisLocation = 'origin';
   ax.YAxisLocation = 'origin';
      
   end
    
   
elseif(k==2)
%% Bairstow Method   
    
    degree=input('Enter the degree of polynomial:');
    coff=input('Enter the coeff in Increasing Order in []:');
    r=input('Enter First Guess :');
    s=input('Enter Second Guess :');
    maxit=input('Enter Maximum no.of Iterations :');
    errmax=input('Enter the Maximum Relative approximate error in(%) :');
    
	scoff=coff;
	sdegree=degree;
    arrc=[];
    d=degree;
    arrb=[];
    flag=0;
    count=0;
    
    while(degree>0)
    
    for j=1:maxit
	
    d=degree;
    arrb(d+1)=coff(d+1);
    arrb(d)=coff(d)+r*arrb(d+1);
    d=d-2;
    
    while(d>=0)
          arrb(d+1)=coff(d+1) + r*arrb(d+2) + s*arrb(d+3);
          d=d-1; 
    end
    
    d=degree;
    arrc(d+1)=arrb(d+1);
    arrc(d)=arrb(d)+ r*arrc(d+1);    
    d=d-2;
    
     while(d>=0)
          arrc(d+1)=arrb(d+1) + r*arrc(d+2) + s*arrc(d+3);
          d=d-1; 
     end
     
     ds = ( arrb(1)*arrc(3) - arrb(2)*arrc(2) )/( arrc(4)*arrc(2) - arrc(3)*arrc(3) );
     dr = ( arrb(1)*arrc(4) - arrb(2)*arrc(3) )/( arrc(3)*arrc(3) - arrc(2)*arrc(4) );
     r=r+dr;
     s=s+ds;
     
     err_r = abs(dr/r)*100 ;
     err_s = abs(ds/s)*100 ;
     
     if ( errmax > err_s ||  errmax > err_r )
         
         root1 = ( r + sqrt(r*r + 4*s) )/2;
         root2 = ( r - sqrt(r*r + 4*s) )/2;
         fprintf('Solution is ');
         
         if(count==0)
             fprintf(' %f %f ',root1,root2);
         end
         count=1;
         flag=1;
         break; 
     end
	 
    end
	
    if flag==0
        break;
    end
	
     for i=1:degree-1
         coff(i)=arrb(i+2);
     end
	 
     degree=degree-2;
	 
     if degree == 2
         root1 = ( -(arrb(4)) + sqrt(arrb(4)*arrb(4) - 4*arrb(3)*arrb(5)) )/(2*arrb(5));
         root2 = ( -(arrb(4)) - sqrt( arrb(4)*arrb(4) - 4*arrb(3)*arrb(5) ))/(2*arrb(5));
         fprintf('%f %f\n ',root1,root2);
         break;
		 
     elseif degree == 1
         root = -arrb(3)/arrb(4);
         fprintf('%f\n ',root);
         break;
		 
     end      
         
    end 
	
	f=zeros(1,101);
   
	
    for k=-50:50
	
        for i=1:sdegree+1
	
            f(k+51)=f(k+51)+scoff(i)*((k)^(i-1));
	 
        end
	   
    end
	   
	 x=-50:50;  
	plot(x,f);grid on;xlabel('x');ylabel('f(x)');title('f(x)vs x');
    ax=gca;
   ax.XAxisLocation = 'origin';
   ax.YAxisLocation = 'origin';
    
%%
else

    fprintf('Enter Valid Number!\n');
    
end