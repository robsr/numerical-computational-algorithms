str = input('Enter the Required Function:','s');
f=inline(str);

j=input('Press 1 for Bisection Method\nPress 2 for False Position Method\nPress 3 for Modified False Position Method\nPress 4 for Newton Raphson Method\nPress 5 for Secant Method\n');
%% Bisection Method
if (j == 1)
    
   xl=input('Enter the lower Guess Value :');xu=input('Enter the upper Guess Value :');
   maxit=input('Enter Maximum no.of Iterations :');errmax=input('Enter the Maximum Relative approximate error in(%) :');
   
   err=[];
   flag=0;
   show=0;
   if (f(xu)*f(xl)>0)
           fprintf('No Root exists between this interval\n');
   elseif (f(xl)==0)
           fprintf('Root of given Equation is :%f\n',xl);
           show=1;
           flag=1;
   elseif(f(xu)==0)
           fprintf('Root of given Equation is :%f\n',xu);
           show=1;
           flag=1;
   else
        xr=(xl+xu)/2;
        pre=xr;
        for i=1:maxit
            if(f(xr) == 0)
                fprintf('Root of given Equation is :%f\n\nNo Error in the Solution\n',xr);
                flag=1;
                show=1;
                break;
            elseif (f(xl)*f(xr) < 0)
               xu=xr;
           
            else   
               xl=xr;
           end
          xr = (xl+xu)/2;
          curr=xr;
          err(i)=abs((curr-pre)/curr)*100;
          pre=curr;
          
          if(err(i)<errmax)
              fprintf('\nRoot of given Equation is :%f\n',curr);
              flag=1;
              break;
          end
       end
   end
   if flag==0
       fprintf('Solution not found\n');   
   elseif(show==0)    
        it=2:i+1;
        plot(it,err,'-r*');grid on;xlabel('No. of Iterations');ylabel('Relative Approximate Error (in%)');title('Error vs No. of Iterations');
   end
   x=-10:10;
   figure;
   plot(x,f(x),'-b');grid on;xlabel('x');ylabel('f(x)');title('f(x) vs x');
   legend('f(x)');    
   ax=gca;
   ax.XAxisLocation = 'origin';
   ax.YAxisLocation = 'origin';

%% False Position Method

elseif (j==2)
    xl=input('Enter the lower Guess Value :');xu=input('Enter the upper Guess Value :');
   maxit=input('Enter Maximum no.of Iterations :');errmax=input('Enter the Maximum Relative approximate error in(%) :');
   
   err=[];
   flag=0;
   show=0;
   if (f(xu)*f(xl)>0)
           fprintf('No Root exists between this interval\n');
   elseif (f(xl)==0)
           fprintf('Root of given Equation is :%f\n',xl);
           show=1;
           flag=1;
   elseif(f(xu)==0)
           fprintf('Root of given Equation is :%f\n',xu);
           show=1;
           flag=1;
   else
        xr=xu-(f(xu)*(xl-xu))/(f(xl)-f(xu));
        pre=xr;
        for i=1:maxit
            if(f(xr) == 0)
                fprintf('Root of given Equation is :%f\n\nNo Error in the Solution\n',xr);
                show=1;
                flag=1;
                break;
            elseif (f(xl)*f(xr) < 0)
               xu=xr;
           else
               xl=xr;
           end
          xr=xu-(f(xu)*(xl-xu))/(f(xl)-f(xu));
          curr=xr;
          err(i)=abs((curr-pre)/curr)*100;
          pre=curr;
          
          if(err(i)<errmax)
              fprintf('\nRoot of given Equation is :%f\n',curr);
              flag=1;
              break;
          end
       end
   end
   if flag==0
       fprintf('Solution not found\nNo. of maximum iterations provided are maybe less\nTry with more no. of iterations\n');
   elseif(show==0)
        it=2:i+1;
        plot(it,err,'-r*');xlabel('No. of Iterations');grid on;ylabel('Relative Approximate Error (in%)');title('Error vs No. of Iterations');
   end
   x=-10:10;
   figure;
   plot(x,f(x),'-b');grid on;xlabel('x');ylabel('f(x)');title('f(x) vs x');
   legend('f(x)');    
   ax=gca;
   ax.XAxisLocation = 'origin';
   ax.YAxisLocation = 'origin';
   
%% Modified False Position Method      

elseif (j==3)
    
     xl=input('Enter the lower Guess Value :');xu=input('Enter the upper Guess Value :');
   maxit=input('Enter Maximum no.of Iterations :');errmax=input('Enter the Maximum Relative approximate error in(%) :');
   
   err=[];
   show=0;
   flag=0;
   if (f(xu)*f(xl)>0)
           fprintf('No Root exists between this interval\n');
   elseif (f(xl)==0)
           fprintf('Root of given Equation is :%f\n',xl);
           flag=1;
           show=1;
   elseif(f(xu)==0)
           fprintf('Root of given Equation is :%f\n',xu);
           flag=1;
           show=1;
   else
        xr=xu-(f(xu)*(xl-xu))/(f(xl)-f(xu));
        pre=xr;
        for i=1:maxit
            if(f(xr) == 0)
                fprintf('Root of given Equation is :%f\n\nNo Error in the Solution\n',xr);
                flag=1;
                show=1;
                break;
            elseif (f(xl)*f(xr) < 0)
               prexu=xu;
               xu=xr;
            else
               prexl=xl;
               xl=xr;
            end
            
          xr=xu-( (f(xu)/2)*(xl-xu) )/(f(xl)-(f(xu)/2));  
          curr=xr;
          err(i)=abs((curr-pre)/curr)*100;
          pre=curr;
          
          if(err(i)<errmax)
              fprintf('\nRoot of given Equation is :%f\n',curr);
              flag=1;
              break;
          end
       end
   end
   if flag==0
       fprintf('Solution not found\nNo. of maximum iterations provided are maybe less\nTry with more no. of iterations\n');
   elseif(show==0)
        it=2:i+1;
        plot(it,err,'-r*');xlabel('No. of Iterations');grid on;ylabel('Relative Approximate Error (in%)');title('Error vs No. of Iterations');
   end
   x=-10:10;
   figure;
   plot(x,f(x),'-b');grid on;xlabel('x');ylabel('f(x)');title('f(x) vs x');
   legend('f(x)');    
   ax=gca;
   ax.XAxisLocation = 'origin';
   ax.YAxisLocation = 'origin';
   
    
%% Newton-Raphson Method
% Derivative is Calculated using Central Difference Method

elseif (j==4)
    
    xg=input('Enter the Guess Value :');
   maxit=input('Enter Maximum no.of Iterations :');errmax=input('Enter the Maximum Relative approximate error in(%) :');
   
   err=[];
   show=0;
   flag=0;
   if (f(xg)==0)
           fprintf('Root of given Equation is :%f\n',xg);
           show=1;
           flag=1;
   else
        xr=xg-f(xg)*.000002/(f(xg+.000001)-f(xg-.000001));
        pre=xr;
        for i=1:maxit
            if (f(xr)==0)
                fprintf('Root of given Equation is :%f\n\nNo Error in the Solution\n',xr);
                flag=1;
                show=1;
                break;
            else
                xr = xr-f(xr)*.000002/(f(xr+.000001)-f(xr-.000001)) ;
                curr=xr;
                err(i)=abs((curr-pre)/curr)*100;
                pre=curr;
            end
          if(err(i)<errmax)
              fprintf('\nRoot of given Equation is :%f\n',curr);
              flag=1;
              break;
          end
       end
   end
   if flag==0
       fprintf('Solution not found\nNo. of maximum iterations provided are maybe less\nTry with more no. of iterations\n');
   elseif(show==0)
        it=2:i+1;
        plot(it,err,'-r*');xlabel('No. of Iterations');grid on;ylabel('Relative Approximate Error (in%)');title('Error vs No. of Iterations');
   end
   x=-10:10;
   figure;
   plot(x,f(x),'-b');grid on;xlabel('x');ylabel('f(x)');title('f(x) vs x');
   legend('f(x)');    
   ax=gca;
   ax.XAxisLocation = 'origin';
   ax.YAxisLocation = 'origin';
   
     
%% Secant Method      
elseif(j==5)
    
    xl=input('Enter the lower Guess Value :');xu=input('Enter the upper Guess Value :');
   maxit=input('Enter Maximum no.of Iterations :');errmax=input('Enter the Maximum Relative approximate error in(%) :');
   
   err=[];
   show=0;
   flag=0;
   if (f(xl)==0)
           fprintf('Root of given Equation is :%f\n',xl);
           show=1;
           flag=1;
   elseif(f(xu)==0)
           fprintf('Root of given Equation is :%f\n',xu);
           show=1;
           flag=1;
   else
        xr = xl-(f(xl)*(xl-xu))/(f(xl)-f(xu));
        pre=xr;
        for i=1:maxit
            if (f(xr)==0)
                fprintf('Root of given Equation is :%f\n\nNo Error in the Solution\n',xr);
                show=1;
                flag=1;
                break;
            else
                xr = xr-(f(xr)*(xr-xl))/(f(xr)-f(xl));
                curr=xr;
                err(i)=abs((curr-pre)/curr)*100;
                pre=curr;
            end
          if(err(i)<errmax)
              fprintf('\nRoot of given Equation is :%f\n',curr);
              flag=1;
              break;
          end
       end
   end
         
   if flag==0
       fprintf('Solution not found\nNo. of maximum iterations provided are maybe less\nTry with more no. of iterations\n');
   elseif(show==0)
        it=2:i+1;
        plot(it,err,'-r*');xlabel('No. of Iterations');grid on;ylabel('Relative Approximate Error (in%)');title('Error vs No. of Iterations');
   end
   x=-10:10;
   figure;
   plot(x,f(x),'-b');grid on;xlabel('x');ylabel('f(x)');title('f(x) vs x');
   legend('f(x)');    
   ax=gca;
   ax.XAxisLocation = 'origin';
   ax.YAxisLocation = 'origin';
          
%%    
else
    fprintf('Enter a valid no.\n');
end
      