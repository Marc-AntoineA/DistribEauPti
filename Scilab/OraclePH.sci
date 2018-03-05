
 
 function [F, G, H, ind] = OraclePH(qc, ind)
     
     F = 0
     G = 0
     H = 0
     
     if ind==2 then
        //-----------valeur du critère au point qc
        //produit scalaire 1
        terme1 = q0 + B*qc
        //produit scalaire 1
        terme2 = r.*(q0 + B*qc).*abs(q0 + B*qc)
        F = (1/3)*terme1'*terme2
        //produit scalaire 2
        terme1 = pr
        terme2 = Ar*(q0 + B*qc)
        F = F + terme1'*terme2
    end
    
    if ind==3 then
        //-----------vecteur des dérivées du critère par rapport à qc
        q = q0 + B*qc
        G = B'*(r.*q.*abs(q) + Ar'*pr)
    end
    
    if ind==4 then
        //-----------valeur du critère au point qc
        //produit scalaire 1
        terme1 = q0 + B*qc
        //produit scalaire 1
        terme2 = r.*(q0 + B*qc).*abs(q0 + B*qc)
        F = (1/3)*terme1'*terme2
        //produit scalaire 2
        terme1 = pr
        terme2 = Ar*(q0 + B*qc)
        F = F + terme1'*terme2
        
        //-----------vecteur des dérivées du critère par rapport à qc
        q = q0 + B*qc
        G = B'*(r.*q.*abs(q) + Ar'*pr)
    end
    
    if ind==5 then
        //matrice des dérivées secondes
        q = q0 + B*qc
        H = B'*diag(r.*(q + abs(q)))*B
    end
    
    if ind==6 then
        //-----------vecteur des dérivées du critère par rapport à qc
        q = q0 + B*qc
        G = B'*(r.*q.*abs(q) + Ar'*pr)
        
        //matrice des dérivées secondes
        H = B'*diag(r.*(q + abs(q)))*B
    end
    
    if ind==7 then
        //-----------valeur du critère au point qc
        //produit scalaire 1
        terme1 = q0 + B*qc
        //produit scalaire 1
        terme2 = r.*(q0 + B*qc).*abs(q0 + B*qc)
        F = (1/3)*terme1'*terme2
        //produit scalaire 2
        terme1 = pr
        terme2 = Ar*(q0 + B*qc)
        F = F + terme1'*terme2
        
        //-----------vecteur des dérivées du critère par rapport à qc
        q = q0 + B*qc
        G = B'*(r.*q.*abs(q) + Ar'*pr)
        
        //matrice des dérivées secondes
        H = B'*diag(r.*(q + abs(q)))*B
    end
