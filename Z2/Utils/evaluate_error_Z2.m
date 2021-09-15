function [mean_err] = evaluate_error_Z2(z_est, zij_orig, Ind)
        Ind_i = Ind(:,1);
        Ind_j = Ind(:,2);
        m = length(Ind_i);
        zij_est = zeros(1,m);
        for k = 1:m
            i=Ind_i(k); j=Ind_j(k); 
            zij_est(k)=z_est(i)*z_est(j);
        end
        mean_err = mean((1-zij_est.*zij_orig)*0.5);

end
