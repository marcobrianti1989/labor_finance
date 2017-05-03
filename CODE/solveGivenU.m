function [TP,Lp_star,w_star_v,EU_vect,...
        E,V,FirmObj,EnteringW0,EnteringLam_Idx,theta_star] = solveGivenU(...
        CV_tol,Niter,nPhi,nG,sep_pol,sigma,pi_Phi,...
        Phi_grid,BETA,gamma_vect,w_star_pre,U,pi_z,r,K,D,tau,w_star_pre_cons,psi,nZ,init_Prod,b,rra)
  
  %Solve the SPP
  [TP,Lp_star,w_star_v,EU_vect] = solvePareto(CV_tol,Niter,nPhi,nG,sep_pol,sigma,pi_Phi,...
    Phi_grid,BETA,gamma_vect,w_star_pre,U,pi_z,r,K,D,tau,w_star_pre_cons);
  
  %Converting the solution from Lagrange multiplier space to promised value space
  [E,V] = calcEV(TP,nPhi,nG,gamma_vect,U,Lp_star,sep_pol,psi,K,D);
  
  %Solve the search problem
  [FirmObj,EnteringW0,EnteringLam_Idx,theta_star] = ...
    solveSearch(nZ,init_Prod,sigma,E,V,U,BETA,EU_vect,b,rra,gamma_vect,sep_pol);
end