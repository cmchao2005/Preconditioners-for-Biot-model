function r=True_solution(uvksip,dx,dy,x,y,t)
%jv 20191031
global mu lambda alpha
if uvksip==1
    if dx==1
        r= cos(x)*exp(-t);
    elseif dy==1
        r=0;
    else
        r=sin(x)*exp(-t);
    end
elseif uvksip==2
     if dx==1
         r = 0;
    elseif dy==1
         r = cos(y)*exp(-t);
     else
         r = sin(y)*exp(-t);
    end
elseif uvksip==3
     if dx==1
        r=alpha*cos(x+y)*exp(-t)+lambda*sin(x)*exp(-t);
    elseif dy==1
        r=alpha*cos(x+y)*exp(-t)+lambda*sin(y)*exp(-t);
     else
        r=alpha*sin(x+y)*exp(-t)-lambda*(cos(x)+cos(y))*exp(-t);
     end
else
    if dx==1
        r=cos(x+y)*exp(-t);
    elseif dy==1
        r=cos(x+y)*exp(-t);
     else
        r=sin(x+y)*exp(-t);
     end
    
end
