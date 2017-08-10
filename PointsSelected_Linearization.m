% Calculate the coefficients for linearized points for GM4th Model

function Coefficients = PointsSelected_Linearization(x, RingTrack_length)

alpha = 0.7;
Vehicle_Num = 20;
Coefficients = zeros(Vehicle_Num,4);

    Coefficients(1,1) = -alpha*x(1+Vehicle_Num)*(x(2*Vehicle_Num)-x(1+Vehicle_Num))/(RingTrack_length-x(1)+x(Vehicle_Num));
    Coefficients(1,2) = alpha*x(1+Vehicle_Num)*(x(2*Vehicle_Num)-x(1+Vehicle_Num))/(RingTrack_length-x(1)+x(Vehicle_Num));
    Coefficients(1,3) = alpha*x(1+Vehicle_Num)/(RingTrack_length-x(1)+x(Vehicle_Num));
    Coefficients(1,4) = alpha*(x(1-1+Vehicle_Num)-2*x(1+Vehicle_Num))/(RingTrack_length-x(1)+x(Vehicle_Num));


for i = 2:Vehicle_Num
    
    Coefficients(i,1) = - alpha*x(i+Vehicle_Num)*(x(i-1+Vehicle_Num)-x(i+Vehicle_Num))/(x(i-1)-x(i));
    Coefficients(i,2) = alpha*x(i+Vehicle_Num)*(x(i-1+Vehicle_Num)-x(i+Vehicle_Num))/(x(i-1)-x(i));
    Coefficients(i,3) = alpha*x(i+Vehicle_Num)/(x(i-1)-x(i));
    Coefficients(i,4) = alpha*(x(i-1+Vehicle_Num)-2*x(i+Vehicle_Num))/(x(i-1)-x(i));
end


