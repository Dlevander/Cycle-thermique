function [X,h_90] = SoutirageTEST(h6,h7,h_80,h9,h_100,nsout,d)
% resolution iterative du systeme de n equation a n inconnues
h_90_new = 2000; % valeur initiale de h_90
err = 10000;
while err >=1e-4
    h_90 = h_90_new;
    A = zeros(nsout,nsout);
    B = zeros(nsout,1);
    
    A(1,1:nsout) = h9(1)-h_90; % ligne 1
    B(1) = -(h9(1)-h_90);
    A(d,1:d-1) = h7(d)-h9(d-1); % ligne d
    B(d) = -h7(d)+h9(d-1);
    
    for i=1:nsout % i ligne
        for j=1:nsout % j colonne
            if j==i % sur la diagonale
                A(j,j) = A(j,j) + (h7(j)-h6(j));
            end
            if j>i
                A(i,j)= A(i,j) + (h7(i)-h7(i+1));
            end
            if i ~=d && i~=1
                A(i,j) = A(i,j)+(h9(i)-h9(i-1));
            end
        end
        if i~=d && i~=1
            B(i) = -(h9(i)-h9(i-1));
        end
    end
    A(1:d-1,d:nsout)= 0; 
    % Resolution
    X = A\B; %fractions soutiree
    
    %%%%%%%%% Calcul etat 90 %%%%%%%%%%
    h_90_new = h_80 + (h_100 - h7(1)) * sum( X(1:d-1) )/( 1+ sum( X(1:d-1) ));
    err = abs(h_90-h_90_new);
end
end