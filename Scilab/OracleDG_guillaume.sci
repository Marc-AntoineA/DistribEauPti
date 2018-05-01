function [y] = g(x)
    y = x./sqrt(abs(x))
endfunction

function [F, G, ind] = OracleDG_guillaume(lambda, ind)
    z = (Ar'*pr + Ad'*lambda)./r;
    qdieze = -g(z)
    F=0;
    G=0;
    if ind==2 then
        F=-1/3*qdieze'*(r.*qdieze.*abs(qdieze))-pr'*Ar*qdieze-lambda'*(Ad*qdieze-fd)
        G=0
    end
    if ind==3 then
        F=0
        G=-Ad*qdieze+fd
    end
    if ind==4 then
        F=-1/3*qdieze'*(r.*qdieze.*abs(qdieze))-pr'*Ar*qdieze-lambda'*(Ad*qdieze-fd)
        G=-Ad*qdieze+fd
    end
    
    //disp(F);
    //disp(G);
    //disp(qdieze);
    disp(norm(z));
    disp(" ");
endfunction
