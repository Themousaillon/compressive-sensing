
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
        zmk=pinv(AS)*x
        [B,h]=gsort(zmk,"g","d")
        a(h(1:s)) = zmk(h(1:s))
        R = x-D*a
    end
    y = a
endfunction



