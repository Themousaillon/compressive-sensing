clear 

// chargement des données
data = csvRead("/Users/nbreizh/Documents/compressive sensing/tp/data.csv")

[M,N] = size(data)

// normalization des données
for i=1:N
    data(:,i) = data(:,i)/norm(data(:,i))
end

function y = calculCoeff(alpha, p, epsilon)
    n = length(alpha)
    
    W = zeros(n, n)
    for i=1:n
        W(i,i) = (abs(alpha(i))^2 + epsilon)^(p/2 -1)
    end
    y = W
endfunction


function y = irls(X, D, p, itermax, epsilon)
    i = 0
    alphak1 = D'*pinv(D*D')*X
    W = calculCoeff(alphak1, p,  epsilon)
    Q=inv(W'*W)
    alphak = Q*D'*pinv(D*Q*D')*X
    while (i < itermax)
        if (abs(norm(alphak) - norm(alphak1)) > sqrt(epsilon)/100) then
            alphak1 = alphak
            W = calculCoeff(alphak, p, epsilon)
            Q=inv(W'*W)
            alphak=Q*D'*pinv(D*Q*D')*X
            i=i+1
         else
             if epsilon > 1e-8 then
                alphak1 = alphak
                epsilon = epsilon/10
                W = calculCoeff(alphak, p, epsilon)
                Q=inv(W'*W)
                alphak= Q*D'*pinv(D*Q*D')*X
            end
            i=i+1
         end
     end
     y = alphak
 endfunction 

// p est la norme à utiliser

function [D,Alpha]=KSVD(X,K,L,p,EPS)
    [N,l]=size(X);
    MAX_ITR=round(K/10);
    D=X(:,1:K);
    s=sqrt(diag(D'*D));
    for i=1:K
        D(:,i)=D(:,i)/s(i);
    end
    Alpha=zeros(K,l); 
    for j=1:L
        for i_vect=1:l
            Alpha(:,i_vect)=irls(X(:,i_vect), D, p, MAX_ITR, EPS);
         end
        for i_col=1:K
            idx_k=find(Alpha(i_col,:)<>0);
            if length(idx_k)>0 then l
                E_k=X-D*Alpha+D(:,i_col)*Alpha(i_col,:);
                Omega=zeros(l,length(idx_k));
                for inz=1:length(idx_k)
                    Omega(idx_k(inz),inz)=1;
                end
                E_kR=E_k*Omega;
                [U,delta,V]=svd(E_kR);
                D(:,i_col)=U(:,1);
                Alpha(i_col,idx_k)=delta(1,1)*V(:,1)';
            else
                g=grand(1,1,"uin",1,l);
               D(:,i_col)=X(:,g)/norm(X(:,g));
            end
        end
    end
endfunction

// Apprentissage à l'aide du stomp

N = 99
K = 100
L = 10
eps=1e-1;
[D,Alpha]=KSVD(data, K,L,0.1, eps);

// On arrondi Alpha pour mieux voir
