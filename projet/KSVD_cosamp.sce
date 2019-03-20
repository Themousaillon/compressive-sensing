clear 

// chargement des données
data = csvRead("/Users/nbreizh/Documents/compressive sensing/tp/data.csv")

[M,N] = size(data)

// normalization des données
for i=1:N
    data(:,i) = data(:,i)/norm(data(:,i))
end


function cj=contribution(R,dj)
    cj=abs(dj'*R)/norm(dj)
endfunction

function y=contriblist(R,D,s)
    [N,M]=size(D)
    CJ=[]
    for i=1:M
        CJ=[CJ,contribution(R,D(:,i))]
    end
    indices=1:length(CJ)
    [B,k]=gsort(CJ,"g","d")
    y=k(1:2*s)
endfunction

function y=cosamp(D,K,x,s)
    [N,M]=size(D)
    a=zeros(M,1)
    R=x
    supp=[]
    for i=1:K
        supp=union(supp, contriblist(R,D,s))
        AS=D(:,supp)
        zmk=AS'*pinv(AS*AS')*x
        [B,h]=gsort(zmk,"g","d")
        a(h(1:s)) = zmk(h(1:s))
        R = x-D*a
    end
    y = a
endfunction



// ss est l'odre de parcimonie souhaitée

function [D,Alpha]=KSVD(X,K,L,ss)
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
            Alpha(:,i_vect)=cosamp(D, MAX_ITR, X(:,i_vect), ss);
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
[D,Alpha]=KSVD(data, K,L,50);

// On arrondi Alpha pour mieux voir
