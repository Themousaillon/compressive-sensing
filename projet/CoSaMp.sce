
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




// tests
// exemple 1

D1 = [1/2*sqrt(2) 1/3*sqrt(3) 1/3*sqrt(6) 2/3 -1/3
-1/2*sqrt(2) -1/3*sqrt(3) -1/6*sqrt(6) 2/3 -2/3
0 -1/3*sqrt(3) 1/6*sqrt(6) 1/3 2/3
]

X1 = [4/3-1/2*sqrt(2); 4/3 + 1/2*sqrt(2);2/3]


// exemple 2

D2 = [
1 1 2 5 0 0 3 -2 1 2 2 2
0 -1 -1 1 0 0 5 0 2 2 7 -1
1 1 1 5 1 2 2 1 1 1 1 5
1 5 2 2 5 0 -4 5 1 5 0 0
0 2 2 1 1 0 0 0 0 4 -1 -2
-1 2 2 2 -2 -3 -4 1 1 1 1 0
]



X2 = [-10;-10;1;21;0;9]
//

disp("solution exemple 1 : ")
disp("cosamp : ")
disp(cosamp(D1, 10, X1, 1))


disp("cosamp : ")
disp(D2*cosamp(D2, 10, X2, 1))
