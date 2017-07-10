function Eprof = calcExpectedProf(wStar,aalpha,phi_d_fun,prod_func_type)
    %Calculate expected profits in a given period given the wStar.
    globalDeclaration;
    
    phi_dw  = phi_d_fun(wStar);
    
    Eprof   = 0;
    for ip  = 1:length(phi_vec)
        phi = phi_vec(ip);
        if phi < phi_dw
            prof_point = 0;
        else if phi > phi_dw
                prof_point = (1-tau)*(prodFn(R,phi,aalpha,r,prod_func_type) - wStar);
            end
        end
        Eprof = Eprof + prof_point;
    end
    Eprof = Eprof/length(phi_vec);
end