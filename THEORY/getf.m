function f_discrete = getf(wStar,phi_db,phi_e,phi_d_fun,aalpha,phi_lim,U_given);
    globalDeclaration
    period = 2;
    E2_discrete = calcExpectedUtil(period,wStar,phi_db,phi_e,phi_d_fun,aalpha,phi_lim);
    period = 3;
    E3_discrete = calcExpectedUtil(period,wStar,phi_db,phi_e,phi_d_fun,aalpha,phi_lim);
    den_discrete = 0;
    switch whichCommitment
        case{'perfect'}
            for ip = 1:length(phi_vec)
                phi           = phi_vec(ip);
                if phi < phi_e
                    den_point = (BETA*(E2_discrete + BETA*utilFunc(b,ssigma,1) ...
                            - (1+BETA)*utilFunc(b,ssigma,1)));
                else if phi > phi_e
                        den_point = (BETA*(E2_discrete + BETA*E3_discrete ...
                                - (1+BETA)*utilFunc(b,ssigma,1)));
                    end
                end
                den_discrete = den_discrete + den_point;
            end
            den_discrete = den_discrete/length(phi_vec);
            f_discrete = (U_given - (1+BETA+BETA^2)*utilFunc(b,ssigma,1))/den_discrete;
        case{'limited'}
            for ip = 1:length(phi_vec)
                phi           = phi_vec(ip);
                if phi < phi_lim
                    den_point = (BETA*(E2_discrete + BETA*utilFunc(b,ssigma,1) ...
                        - (1+BETA)*utilFunc(b,ssigma,1)));
                else if phi >= phi_lim
                        den_point = (BETA*(E2_discrete + BETA*E3_discrete ...
                            - (1+BETA)*utilFunc(b,ssigma,1)));
                    end
                end
                den_discrete = den_discrete + den_point;
            end
            den_discrete = den_discrete/length(phi_vec);
            f_discrete = (U_given - (1+BETA+BETA^2)*utilFunc(b,ssigma,1))/den_discrete;
    end
end