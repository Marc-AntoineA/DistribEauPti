
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
    
     titre = "Parametres du gradient a pas fixe";
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
        else
            delta_x = alpha*D;
            delta_G = G - G_precedent;
            y = 1/(delta_G'*delta_x);
            W = (eye(n-md,n-md) - y*delta_x*delta_G') * W_precedent * (eye(n-md,n-md) - y*delta_x*delta_G') + y*delta_x*delta_x';
        end
        
        //    - calcul de la longueur du pas de gradient
        
        [alpha,ok] = Wolke_Skel(alpha,x,D,Oracle);
        if ok==2 then
            alpha = 1;
        end
        
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
    disp('Fin de la methode de gradient a pas fixe')
    disp(cvge)
    
    // - visualisation de la convergence
    
    Visualg(logG,logP,Cout);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    [F,G,H,ind] = Oracle(x,6);
    [F_precedent,G_precedent,H_precedent,ind] = Oracle(x_precedent,3);
    delta_x = x - x_precedent;
    delta_G = G - G_precedent;
    
    I = eye(n-md,n-md);
    y = 1/delta_G'*delta_x;
    W = (I - y*(delta_x*delta_G')) * W_precedent * (y * delta_x*delta_G');
    
    d = -W*G;
    
endfunction
