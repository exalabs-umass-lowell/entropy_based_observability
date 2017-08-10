
% By using renormalization group theory to show the variation under
% different scales for a random Lattice


% Original lattice, being white grid cell for 1, and black grid cell for 0
x = 1:64;
y = 1:64;
[X,Y] = meshgrid(x,y);
%Note that size(Z) is the same as size(x) and size(y)
Z = round(0.8*rand(64));   
pcolor(x,y,Z);
colormap(gray(2)); % make sure its color shows white and black


%% Variation under different scales

N = 3;   % the exponent of maximum ratio of scale over the origin

for i = 1:N
    
    Scale = 2^i;
    
    Len = length(x)/Scale;      % length under specific scale
    Wid = length(y)/Scale;      % width under specific scale
    Pole = zeros(Len, Wid);     % Value at the specific point
    for j = 1:Len
        for k = 1:Wid
            %  Here calculate the mean value of the cell area in the
            %  original graph, and take 0.2 as the criteria to determine
            %  it's white or black
            
            Epsilon = 0.2;
            if mean(Z(Scale*(j-1)+1:Scale*j,Scale*(k-1)+1:Scale*(k))) >=Epsilon   
                Pole(j,k) = 1;
            else
                Pole(j,k) = 0;
       
            end
        end
    end
    
figure(i+1);
pcolor(1:Len,1:Wid,Pole);
colormap(gray(2))
end



