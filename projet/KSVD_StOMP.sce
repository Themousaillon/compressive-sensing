
clear 

// chargement des données
data = csvRead("/Users/nbreizh/Documents/compressive sensing/tp/data.csv")

[N,M] = size(data)

// normalization des données
for i=1:M
    data(:,i) = data(:,i)/norm(data(:,i))
end


function y = majIndices(p, alpha)
    y = []
    for i=1:length(p)
        if alpha(p(i)) <> 0 then
            y = [y, p(i)]
        end
    end
endfunction

function cj=contribution(R,dj)
    cj=abs(dj'*R)/norm(dj)
endfunction

function y=contriblist(R,D)
    [N,M]=size(D)
    CJ=[]
    for i=1:M
        CJ=[CJ,contribution(R,D(:,i))]
    end
    y= CJ
endfunction

function y= selection_atomes(dic, R, t)
    [M,N] = size(dic)
    S = t*norm(R)/sqrt(N)
    CL = contriblist(R,dic)
    y = find(CL>S)
endfunction

function y = stOMP(dic, X, t, K, epsilon)
    [M,N] = size(dic)
    alpha = zeros(N,1)
    R=X
    i=1
    // on stockera les indices des atomes choisis dans p
    p = []
    while (i<=K) && (norm(dic*alpha-X)>epsilon)
        lambda = selection_atomes(dic, R, t)
        p=union(p,lambda)
        phi = dic(:, p)
        zmK = phi'*pinv(phi*phi')*X
        alpha(p) = zmK
        R=X-dic*alpha
        p = majIndices(p, alpha)
        i = i+1
    end
    y = alpha
endfunction

// t est le critère de seuillage à passer à l'stomp

function [D,Alpha]=KSVD(X,K,L,t,EPS)
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
            Alpha(:,i_vect)=stOMP(D, X(:,i_vect),t, MAX_ITR, EPS);
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
eps=1e-6;
[D,Alpha]=KSVD(data, K,L,2.5,eps);
