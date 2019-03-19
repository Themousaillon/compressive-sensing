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
    alphak1 = D'*inv(D*D')*X
    W = calculCoeff(alphak1, p,  epsilon)
    Q=inv(W'*W)
    alphak = Q*D'*inv(D*Q*D')*X
    while (i < itermax)
        if (abs(norm(alphak) - norm(alphak1)) > sqrt(epsilon)/100) then
            alphak1 = alphak
            W = calculCoeff(alphak, p, epsilon)
            Q=inv(W'*W)
            alphak=Q*D'*inv(D*Q*D')*X
            i=i+1
         else
             if epsilon > 1e-8 then
                alphak1 = alphak
                epsilon = epsilon/10
                W = calculCoeff(alphak, p, epsilon)
                Q=inv(W'*W)
                alphak= Q*D'*inv(D*Q*D')*X
            end
            i=i+1
         end
     disp(epsilon)
     end
     y = alphak
 endfunction 
