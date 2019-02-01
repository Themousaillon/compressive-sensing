clear 


function y= selection_atome(dic, R)
    s = size(dic)
    m=1
    maxi=abs(dic(:,1)'*R)/norm(dic(:,1))
    for i=2:s(2)
        tmp = abs(dic(:,i)'*R)/norm(dic(:,i))
        if tmp > maxi then
            maxi=tmp
            m=i
        end
    end
    y=m
endfunction


function y = ortogonal_MP(dic, X, K, epsilon)
    s = size(dic)
    alpha = zeros(s(2),1)
    R=X
    i=1
    // on stockera les indices des atomes choisis dans p
    p = []
    while (i<=K) && (norm(dic*alpha-X)>epsilon)
        m = selection_atome(dic, R)
        p=[p;m]
        phi(:, p) = dic(:, p)
        zmK = X'*phi(:, p)*inv(phi(:, p)'*phi(:, p))
        alpha(p) = zmK'
        R=X-dic*alpha
        i=i+1
    end
    y = alpha
endfunction

// Vous devez executer au préalable la fonction OMP ou l'appeler : "exec(chemin de la fonction)"
///////////////////////////////////////////////////////////////////////////////////////////////////////////
function [D,Alpha]=KSVD(X,K,L,EPS)
    //En entrée on a une matrice X constituée des signaux d'apprentissage, chaque colonne corespond à un vecteur d'apprentissage : On a l signaux d'apprentissage et chaque signal est de dimension N
//K est le nombre de colonnes de dictionnaire qu'on veut apprendre
// On doit avoir l plus grand que K
//L est le nombre de fois de mise à jour de dictionnaire 
    [N,l]=size(X);
    MAX_ITR=round(K/10);// le nombre d'iteration maximal à effectuer pour l'OMP -- pas besoin de faire beaucoup d'itérations
    D=X(:,1:K); // On initialise le dictionnaire en considérant les premières colonnes des veceturs d'apprentissage
    s=sqrt(diag(D'*D));
    for i=1:K
        D(:,i)=D(:,i)/s(i); // On normalise les colonnes du dictionnaire
    end
    Alpha=zeros(K,l); // On cherche l solutions parcimonieuses pour les vecteurs d'apprentissage, les alpha_i sont de dimension égale à K 
    for j=1:L //Répéter les mises à jour L fois
        for i_vect=1:l
            Alpha(:,i_vect)=ortogonal_MP(D, X(:,i_vect),MAX_ITR, EPS); //OMP sur chaque vecteur d'apprentissage et on renvoie la solution parcimonieuse
         end
        //D=D0;
        for i_col=1:K // Effectuer une SVD sur chaque colonne==> K colonnes issues de SVD
            idx_k=find(Alpha(i_col,:)<>0); // Déterminer le support de la ligne i_col de Alpha
            if length(idx_k)>0 then // cas où support est non vide ==> card(supp) non nul
                E_k=X-D*Alpha+D(:,i_col)*Alpha(i_col,:);
                Omega=zeros(l,length(idx_k));
                for inz=1:length(idx_k)
                    Omega(idx_k(inz),inz)=1;
                end
                E_kR=E_k*Omega;
                [U,delta,V]=svd(E_kR);
                D(:,i_col)=U(:,1);
                Alpha(i_col,idx_k)=delta(1,1)*V(:,1)';
            else // si card (supp)=0, on choisit de façon aléatoire une colonne de la matrice X
                g=grand(1,1,"uin",1,l); // choix aléatoire d'un indice de 1 à l
               D(:,i_col)=X(:,g)/norm(X(:,g)); // La colonne du dictionnaire est donc égal à une colonne de X qu'on normalise
            end
        end
    end
endfunction
X=grand(30,75,"uin",1,4);// j'ai généré une matrice X d'apprentissage 55 vecteurs(colonnes) de taille 30
eps=1e-6;
[D,Alpha]=KSVD(X, 30,100,eps);
