function solve_test

x = fsolve(@myfun,[0 1], optimoptions('fsolve','Display','iter'))

disp('check the nuerical solution, make sure the value is close to 0');
myfun(0.719139328218676)

function F = myfun(x)
    F = x.^5-2*x.^3 + 3*x.^2-1;        
end

end
