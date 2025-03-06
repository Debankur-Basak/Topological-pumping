










f=0;

for n=1:2:7
    for m=1:2:7

f=f+sin(n*pi/2)*sin(m*pi/2)*(n*m*pi^2/4)*(exp(-pi*sqrt(n^2+m^2)/2)+((1-exp(-pi*sqrt(n^2+m^2)))*sinh(pi*sqrt(n^2+m^2)/2)/sinh(pi*sqrt(n^2+m^2))));
    end
end

disp(f);