function XMASSFLOW = Soutirage(h6,h7,h_80,h9,h_100,nsout,d)
% resolution iterative du systeme de n equation a n inconnues
h_90_new = 2000; % valeur initiale de h_90
err = 10000;
while err >=0.001
    h_90 = h_90_new;
    A = zeros(nsout,nsout);
    B = zeros(nsout,1);
    %remplissage de la matrice A et B
    %premiere ligne particuliere
    A(1,1:d-1) = A(1,1:d-1) + ( h9(1) - h_90);
    A(1,1) = A(1,1) - ( h7(1) -h6(1) );
    A(1,2:d-1) = A(1,2:d-1) - ( h7(1) - h7(2) );
    
    B(1) = (h9(1)-h_90);
    
    for i=2:nsout % i les lignes de la matrice
        if i < nsout
            deltaH7  = h7(i) - h7(i+1);
        end
        deltaH9  = h9(i) - h9(i-1);
        deltaH76 = h7(i) - h6(i);
        if i < d
            A(i,1:d-1) = A(i,1:d-1) + deltaH9;
            A(i,i) = A(i,i) - deltaH76; % remplissage diagonale
            A(i,i+1:d-1) = A(i,i+1:d-1) - deltaH7;
            B(i) = deltaH9;
        end
        if i > d
            A(i,:) = A(i,:) + deltaH9;
            A(i,i) = A(i,i) - deltaH76; % remplissage diagonale
            if i < nsout
                A(i,i+1:nsout) = A(i,i+1:nsout) - deltaH7;
            end
            B(i) = deltaH9;
        end
    end
    % remplissage ligne d particuliere
    A(d,1:d-1) = A(d,1:d-1) + h7(d) - h9(d-1);
    A(d,d) = h7(d) - h6(d);
    A(d,d+1:nsout) = A(d,d+1:nsout) + h7(d) - h7(d+1);
    %Ax=B
    B(d) =-( h7(d) - h9(d-1));
    % Resolution
    
    X = A\B; %fractions soutiree
    XMASSFLOW = X;
   
    %%%%%%%%% Calcul etat 90 %%%%%%%%%%
    h_90_new = h_80 + (h_100 - h7(1)) * sum( XMASSFLOW(1:d-1) )/( 1+ sum( XMASSFLOW(1:d-1) ));
    err = abs(h_90-h_90_new);
end
    
end