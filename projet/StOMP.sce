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
        zmK = phi'*inv(phi*phi')*X
        alpha(p) = zmK
        R=X-dic*alpha
        p = majIndices(p, alpha)
        i = i+1
    end
    y = alpha
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
disp("stOMP : ")
disp(stOMP(D1, X1, 2.5, 10, 0.001))


disp("stOMP : ")
disp(D2*stOMP(D2, X2, 2.5, 10, 0.0001))
