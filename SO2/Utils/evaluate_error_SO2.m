function [mean_err] = evaluate_error_SO2(theta_est, thetaij_orig, Ind)
        Ind_i = Ind(:,1);
        Ind_j = Ind(:,2);
        m = length(Ind_i);
        thetaij_est = zeros(1,m);
        for k = 1:m
            i=Ind_i(k); j=Ind_j(k); 
            thetaij_est(k)=theta_est(i)-theta_est(j);
        end
        thetaij_est = mod(thetaij_est+2*pi, 2*pi);
        err_vec = mod(thetaij_est-thetaij_orig+2*pi, 2*pi)/pi;
        mean_err = mean(min(err_vec,2-err_vec));

end