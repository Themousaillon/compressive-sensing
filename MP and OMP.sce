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


function y = matching_pursuit(dic, X, K)
    s = size(dic)
    alpha = zeros(s(2),1)
    R=X
    for i=1:K
        m = selection_atome(dic, R)
        zmK = ((dic(:,m)'*R))/(norm(dic(:,m))^2)
        R = R - zmK*dic(:,m)
        alpha(m) = alpha(m) + zmK
    end
    y = alpha
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
disp("MP : ")
disp(matching_pursuit(D1, X1, 10))
disp("OMP : ")
disp(ortogonal_MP(D1, X1, 10, 0.001))


disp("solution exemple 2 : ")
disp("MP : ")
sol = matching_pursuit(D2, X2, 6)
disp(matching_pursuit(D2, X2, 1000))
disp("OMP : ")
disp(ortogonal_MP(D2, X2, 1000, 0.0000001))
solotho = ortogonal_MP(D2, X2, 6, 0.0000001)
