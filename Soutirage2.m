function X = Soutirage2(h6,h7,h_80,h9,h_100,nsout,d)
% resolution iterative du systeme de n equation a n inconnues
h_90_new = 20; % valeur initiale de h_90
err = 10000;
deltaH76 = diag(h7-h6);
deltaH7 = h7(1:nsout-1)-h7(2:nsout);
while err >=0.001
    h_90 = h_90_new;
    A = zeros(nsout,nsout);
    if d ~= 1
        %Remplissage A coin avant ligne d , matrice taille d-1 x d-1
        H9 = [h_90 h9(1:nsout-1)];
        deltaH9 = h9-H9;
        A(1:d-1,1:d-1) = A(1:d-1,1:d-1)+deltaH9(1:d-1)'*ones(1,d-1)-deltaH76(1:d-1,1:d-1);
        if d > 2
            for i=1:d-2
                A(1:d-1,1:d-1) = A(1:d-1,1:d-1) -diag(deltaH7(1:d-1-i),i);
            end
        end
    end
    if d < nsout
        %Remplissage A coin apres ligne d , matrice taille nsout-d x nsout-d
        A(d+1:nsout,d+1:nsout) =  A(d+1:nsout,d+1:nsout)-deltaH76(d+1:nsout,d+1:nsout);
        if d < nsout-2
            for i=1:nsout-d-1
                A(d+1:nsout,d+1:nsout) = A(d+1:nsout,d+1:nsout) -diag(deltaH7(d+1:nsout-i),i);
            end
        end
        %Remplissage apres ligne d
        A(d+1:nsout,:) = A(d+1:nsout,:)+deltaH9(d+1:nsout)'*ones(1,nsout);
    end
    
    %Remplissage ligne d
    A(d,1:nsout)= h7(d);
    A(d,d) = A(d,d)-h6(d); 
    if d < nsout
        A(d,d+1:nsout) = A(d,d)-h7(d+1); % colonne a droite de d
    elseif d > 1
        A(d,1:d-1) = A(d,1:d-1)-H9(d-1); %colonne a gauche de d
    end
    %remplissage B
    B = -deltaH9';
    B(d) = B(d) + h9(d) - h7(d); 
    % Resolution
    X = A\B; %fractions soutiree
    
    %%%%%%%%% Calcul etat 90 %%%%%%%%%%
    h_90_new = h_80 + (h_100 - h7(1)) * sum( X(1:d-1) )/( 1+ sum( X(1:d-1) ));
    err = abs(h_90-h_90_new);
end