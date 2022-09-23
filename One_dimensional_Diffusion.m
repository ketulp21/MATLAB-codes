clear all
close all
clc

%% Defining the mesh 

n_points=5;
dom_size=1;
h=dom_size/(n_points-1);

%% initialising the problem
y(1)=0;
y(5)=1;

y_new(1)=0;
y_new(n_points)=1;

error_mag=1;
error_req=1e-6;
iterations=0;


%% calculations

while error_mag>error_req
    for i=2:(n_points-1)
        y_new(i)=(y(i-1)+y(i+1))/2;
    end
    iterations=iterations+1;
    %Calculation of error magnitude
    error_mag=0;
    for i=2:n_points
        error_mag=error_mag+(abs(y(i)-y_new(i)));
    end
    y=y_new;
end
 %% plotting

 x_dom=((1:n_points)-1)*h;
 figure;
 plot(x_dom,y);