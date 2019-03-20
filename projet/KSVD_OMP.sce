
clear 

// chargement des données
data = csvRead("/Users/nbreizh/Documents/compressive sensing/tp/data.csv")

[M,N] = size(data)

// normalization des données
for i=1:N
    data(:,i) = data(:,i)/norm(data(:,i))
end


function y= selection_atome(dic, R)
    [M,N] = size(dic)
    m=1
    maxi=abs(dic(:,1)'*R)/norm(dic(:,1))
    for i=2:N
        tmp = abs(dic(:,i)'*R)/norm(dic(:,i))

        if tmp > maxi then
            maxi=tmp
            m=i
        end
    end
    y=m
endfunction

function y = ortogonal_MP(dic, X, K, epsilon)
    [M,N] = size(dic)
    alpha = zeros(N,1)
    R=X
    i=1
    // on stockera les indices des atomes choisis dans p
    p = []
    while (i<=K) && (norm(dic*alpha-X)>epsilon)
        m = selection_atome(dic, R)
        p=union(p,m)
        phi = dic(:, p)
        zmK = phi'*pinv(phi*phi')*X
        alpha(p) = zmK
        R=X-dic*alpha
        i=i+1
    end
    y = alpha
endfunction

function [D,Alpha]=KSVD(X,K,L,EPS)
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
            Alpha(:,i_vect)=ortogonal_MP(D, X(:,i_vect), MAX_ITR, EPS);
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
eps=1e-8;
[D,Alpha]=KSVD(data, K,L,eps);

// On arrondi Alpha pour mieux voir



