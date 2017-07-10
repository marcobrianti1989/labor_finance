function g_discrete = getg(wStar,aalpha,phi_d_fun,phi_lim,phi_e,prod_func_type);
    globalDeclaration;
    Eprof = calcExpectedProf(wStar,aalpha,phi_d_fun,prod_func_type);
    g_discrete = 0;
    switch whichCommitment
        case{'perfect'}
            for ip = 1:length(phi_vec)
                phi           = phi_vec(ip);
                if phi < phi_e
                    g_point = BETA*Eprof;
                else if phi > phi_e
                        g_point = BETA*Eprof*(1+BETA);
                    end
                end
                g_discrete = g_discrete + g_point;
            end
            g_discrete = g_discrete/length(phi_vec);
        case{'limited'}
            for ip = 1:length(phi_vec)
                phi           = phi_vec(ip);
                if phi < phi_lim
                    g_point = BETA*Eprof;
                else if phi >= phi_lim
                        g_point = BETA*Eprof*(1+BETA);
                    end
                end
                g_discrete = g_discrete + g_point;
            end
            g_discrete = g_discrete/length(phi_vec);
    end