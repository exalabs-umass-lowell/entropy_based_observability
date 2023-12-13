function [Energy, Count] = N_Choose_K (N, K)


% number of vehicles
Nv = K
% number of empty sites
Ne = N-Nv;

Energy = zeros(Nv, 1);
Count = zeros(Nv, 1)
for i = 1:Nv
    Energy(i) = Nv-1+Ne-1-4*i+2;
    Split_Num = nchoosek(Nv-1, i-1);
    
    Split_Num_prev = 0;
    if i ~= Nv
        Split_Num_prev = nchoosek(Nv-1, i);
    end
   
    if(i<2)
        Insert_Num = nchoosek(Ne+1, i);
        Count(i) = Insert_Num + Split_Num_prev;
    % assumptions:  more vehicles than empty sites
    elseif(Ne-1 < i)
        Insert_Num = nchoosek(Ne-1, i-1)*2;
        SameInEnd = nchoosek(Ne-1, i-2);
        Count(i,1) = Split_Num*Insert_Num+SameInEnd*Split_Num_prev;
    else
        Insert_Num = nchoosek(Ne-1, i)+nchoosek(Ne-1, i-1)*2;
        SameInEnd = nchoosek(Ne-1, i-2);
        Count(i,1) = Split_Num*Insert_Num+SameInEnd*Split_Num_prev;
    end
    % consider the case that left end and right end are both vehicles
    
end

end

