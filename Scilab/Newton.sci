function [fopt,xopt,gopt]=Newton(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de gradient a pas variable                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres pour Newton";
   labels = ["Nombre maximal d''iterations";...
             "Valeur du pas de gradient";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1,"vec",1);
   default = ["5000";"0.0005";"0.000001"];
   [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

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
   alpha = alphai;
   kstar = iter;
   for k = 1:iter

//    - valeur du critere, du gradient et du Hessien

      ind = 7;
      [F,G,H] = Oracle(x,ind);

//    - test de convergence

      if norm(G) <= tol then
         kstar = k;
         break
      end

//    - calcul de la direction de descente

      D = -inv(H)*G;

//    - calcul de la longueur du pas de gradient

      [alpha, ok] = Wolfe(1, x, D, Oracle);
      //alpha = 1
//    - mise a jour des variables

      x = x + (alpha*D);

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
   disp('Fin de la methode de Newton a pas variable')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout,5);

endfunction
