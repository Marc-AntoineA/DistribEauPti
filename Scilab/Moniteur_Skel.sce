///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  MONITEUR D'ENCHAINEMENT POUR LE CALCUL DE L'EQUILIBRE D'UN RESEAU D'EAU  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --------------------------------------
// Dimensionnement de l'espace de travail
// --------------------------------------

   stacksize(10000000);

// ------------------------------------------
// Fonctions fournies dans le cadre du projet
// ------------------------------------------

   // Donnees du problemes

   exec('Probleme_R.sce');
   exec('Structures_R.sce');
   
    // Resolutions
    exec('Gradient_F.sci');
    exec('Gradient_V.sci');
    exec('BFGS.sci');
    exec('Polak_Ribiere.sci');
    exec('Wolfe_Skel.sci');
    exec('Optim_Scilab.sci');
    exec('Newton.sci');


   // Affichage des resultats

   exec('Visualg.sci');
   
   // Verification  des resultats
   
   exec('HydrauliqueP.sci');
   exec('HydrauliqueD.sci');
   exec('Verification.sci');

// ------------------------------------------
// Fonctions a ecrire dans le cadre du projet
// ------------------------------------------

   // ---> Charger les fonctions  associees a l'oracle du probleme,
   //      aux algorithmes d'optimisation et de recherche lineaire.
   //
   // Exemple : la fonction "optim" de Scilab
   //
   exec('OraclePG.sci');
   exec('OraclePH.sci');
   exec('OracleDG.sci');
   exec('OracleDH.sci');
   exec('OracleDG_guillaume.sci');
   
   exec('Optim_Scilab.sci');
   titrgr = "Fonction optim de Scilab sur le probleme primal";

   // -----> A completer...
   // -----> A completer...
   // -----> A completer...

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

   // La dimension (n-md) est celle du probleme primal

   xini = 0.1 + 0* rand(n-md,1);
   lambda_ini = 0.1 * rand(md,1);

// ----------------------------
// Minimisation proprement dite
// ----------------------------

   // Exemple : la fonction "optim" de Scilab
   //
   [fopt,xopt,gopt] = Optim_Scilab(OraclePG,xini);

   
   titrgr = "Gradient à pas fixé";
   //Gradient_F(OraclePG, xini, 1);
   titrgr = "Gradient à pas fixé : problème dual";
   Gradient_F(OracleDG, lambda_ini, 2);
   Gradient_F(OracleDG_guillaume, lambda_ini, 2);
   
   titrgr = "Gradient à pas variable";
   //Gradient_V(OraclePG, xini, 3);
   titrgr = "Gradient à pas variable : problème dual";
   //Gradient_V(OracleDG, lambda_ini, 4);
   
   titrgr = "Polak ribiere";
   //Polak_Ribiere(OraclePG, xini, 5);
   titrgr = "Polak ribiere : problème dual";
   //Polak_Ribiere(OracleDG, lambda_ini, 6);
   
   titrgr = "BFGS";
   //BFGS(OraclePG, xini, 7);
   titrgr = "BFGS : problème dual";
   //BFGS(OracleDG, lambda_ini, 8);
   
   titrgr = "Newton";
   //Newton(OraclePH, xini, 9);
   titrgr = "Newton : problème dual";
   //Newton(OracleDH, lambda_ini, 10);
   

// --------------------------
// Verification des resultats
// --------------------------

   [q,z,f,p] = HydrauliqueP(xopt);

   Verification(q,z,f,p);

//
