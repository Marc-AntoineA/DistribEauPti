
function [fopt,xopt,gopt] = BFGS(Oracle,xini)
    
    ///////////////////////////////////////////////////////////////////////////////
    //                                                                           //
    //         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
    //                                                                           //
    //         Methode de gradient a pas fixe                                    //
    //                                                                           //
    ///////////////////////////////////////////////////////////////////////////////


    // ------------------------
    // Parametres de la methode
    // ------------------------
    
     titre = "Parametres de BFGS";
     labels = ["Nombre maximal d''iterations";...
                 "Valeur du pas de gradient";...
                 "Seuil de convergence sur ||G||"];
     typ = list("vec",1,"vec",1,"vec",1);
     default = ["5000";"0.0005";"0.000001"];
     [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);
     alpha = alphai;
    
    // ----------------------------
    // Initialisation des variables
    // ----------------------------
       
     logG = [];
     logP = [];
     Cout = [];
    
     timer();
    
    // -------------------------
    // Boucle sur les iterations
    // -------------------------
    
     x = xini;
    
     kstar = iter;
     for k = 1:iter
    
        //    - valeur du critere et du gradient
        ind = 4;
        [F,G] = Oracle(x,ind);
        //    - test de convergence
        
        if norm(G) <= tol then
          kstar = k;
          break
        end
        
        //    - calcul de la direction de descente
        
        if k==1 then
            D = -G;
            W = eye(n-md,n-md);
        else
            delta_x = alpha*D;
            delta_G = G - G_precedent;
            y = 1/(delta_G'*delta_x);
            A = eye(n-md,n-md) - y*delta_x*delta_G';
            W = A * W_precedent * A + y*delta_x*delta_x';
            D = -W*G;
        end
        //    - calcul de la longueur du pas de gradient
        [alpha,ok] = Wolfe(1,x,D,Oracle);
        
        //    - mise a jour des variables
        
        x = x + (alpha*D);
        G_precedent = G;
        W_precedent = W;
        
        //    - evolution du gradient, du pas et du critere
        
         logG = [ logG ; log10(norm(G)) ];
         logP = [ logP ; log10(alpha) ];
         Cout = [ Cout ; F ];
    
       end
    
    // ---------------------------
    // Resultats de l'optimisation
    // ---------------------------
    
    fopt = F;
    xopt = x;
    gopt = G;
    
    tcpu = timer();
    
    cvge = ['Iteration         : ' string(kstar);...
               'Temps CPU         : ' string(tcpu);...
               'Critere optimal   : ' string(fopt);...
               'Norme du gradient : ' string(norm(gopt))];
    disp('Fin de la methode de BFGS')
    disp(cvge)
    
    // - visualisation de la convergence
    
    Visualg(logG,logP,Cout);
    
endfunction
