function [mean_err] = evaluate_error_SO3(R_est, Rij_orig, Ind)
        Ind_i = Ind(:,1);
        Ind_j = Ind(:,2);
        m = length(Ind_i);
        d = 3;
        Rij_est = zeros(d,d,m);
        for k = 1:m
            i=Ind_i(k); j=Ind_j(k); 
            Rij_est(:,:,k)=R_est(:,:,i)*(R_est(:,:,j)');
        end
        
        MSEVec = zeros(1,m);
        for k = 1:m
            R_tr = trace(Rij_orig(:,:,k)*(Rij_est(:,:,k))');
            MSEVec(k) =  abs(acos((R_tr-1)./2))/pi;
        end
        mean_err = mean(MSEVec);

end