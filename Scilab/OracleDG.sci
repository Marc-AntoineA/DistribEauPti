
function [q_diez] = calc_q_diez(lambda)
    //calcule le q_diez, résultat le la minimisation en q
    z0 = (Ad'*lambda + Ar'*pr) ./ r;
    q_diez = sign(z0) .* sqrt(abs(z0));
    //q_diez = - z ./ sqrt(abs(z));
endfunction

function [L] = lagrange(q, lambda)
    //calcule le lagrangien du problème en q, lambda
    L = q'*(r.*q.*abs(q))/3 + pr'*Ar*q + lambda'*(Ad*q - fd);
endfunction


function [F, G, ind] = OracleDG(lambda, ind)
    
    F = 0
    G = 0
    z0 = (Ad'*lambda + Ar'*pr) ./ r;
    q_diez = - sign(z0) .* sqrt(abs(z0));
    //q_diez = - z0 ./ sqrt(abs(z0));
    
    if ind==2 then
        //-----------valeur de phi en lambda
        [F] = lagrange(q_diez, lambda);
    end
    
    if ind==3 then
        //-----------gradient de phi par rapport à lambda
        G = Ad*q_diez - fd;
    end
    
    if ind==4 then
        //-----------valeur de phi en lambda
        [F] = lagrange(q_diez, lambda);
        
        //-----------gradient de phi par rapport à lambda
        G = Ad*q_diez - fd;
    end
    
    F = -F;
    G = -G;
    
 endfunction
