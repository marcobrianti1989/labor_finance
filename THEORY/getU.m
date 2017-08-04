function U = getU(params,output,phi_cutoff,wageStar,ptheta,aalpha)
  
  E2_discrete   = calcExpectedUtil(params,output,phi_cutoff,wageStar);
  E3_discrete   = getE3(params,aalpha,wageStar);
  limitIntegral = getContinuationProb(phi_cutoff,params.phi_vec);

  U             = params.utilFunc(params.b) ...
    + ptheta*(params.BETA*E2_discrete ...
    + params.BETA^2*(limitIntegral*params.utilFunc(params.b) ...
    + (1-limitIntegral)*E3_discrete)) ...
    + (1-ptheta)*(params.BETA + params.BETA^2)*params.utilFunc(params.b);
end