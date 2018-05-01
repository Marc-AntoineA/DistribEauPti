function [q_diez] = calc_q_diez(lambda)
    //calcule le q_diez, résultat le la minimisation en q
    z = -(Ad'*lambda + Ar'*pr) ./ r;
    q_diez = sign(z) .* z;
endfunction

function L = lagrange(q, lambda)
    //calcule le lagrangien du problème en q, lambda
    L = q'*(r.*q.*abs(q))/3 + pr'*Ar*q + lambda'*(Ad*q - fd);
endfunction


function [F, G, H, ind] = OracleDH(lambda, ind)
     
    F = 0
    G = 0
    H = 0
     
    if ind==2 then
        //-----------valeur de phi en lambda
        q_diez = calc_q_diez(lambda);
        F = lagrange(q_diez, lambda);
    end
    
    if ind==3 then
        //-----------gradient de phi par rapport à lambda
        q_diez = calc_q_diez(lambda);
        G = Ad*q_diez - fd;
    end
    
    if ind==4 then
        //-----------valeur de phi en lambda
        q_diez = calc_q_diez(lambda);
        F = lagrange(q_diez, lambda);
        
        //-----------gradient de phi par rapport à lambda
        G = Ad*q_diez - fd;
    
    end
    
    if ind==5 then
        //----------hessienne de phi par rapport à lambda
        q_diez = calc_q_diez(lambda);
        d = ones(q_diez) ./ (abs(q_diez) .* r);
        H = -Ad*diag(d)*Ad' / 2;
    end
    
    if ind==6 then
        //-----------gradient de phi par rapport à lambda
        q_diez = calc_q_diez(lambda);
        G = Ad*q_diez - fd;
        
        //----------hessienne de phi par rapport à lambda
        d = ones(q_diez) ./ (abs(q_diez) .* r);
        H = -Ad*diag(d)*Ad' / 2;
    end
    
    if ind==7 then
        //-----------valeur de phi en lambda
        q_diez = calc_q_diez(lambda);
        F = lagrange(q_diez, lambda);
        
        //-----------gradient de phi par rapport à lambda
        G = Ad*q_diez - fd;
        
        //----------hessienne de phi par rapport à lambda
        d = ones(q_diez) ./ (abs(q_diez) .* r);
        H = -Ad*diag(d)*Ad' / 2;
    end
    
    F = -F;
    G = -G;
    H = -H;
    
endfunction
