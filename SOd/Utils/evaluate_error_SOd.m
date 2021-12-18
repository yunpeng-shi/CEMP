function [mean_err] = evaluate_error_SOd(R_est, Rij_orig, Ind)

        d = size(R_est,1);
        Ind_i = Ind(:,1);
        Ind_j = Ind(:,2);
        m = length(Ind_i);
        Rij_est = zeros(d,d,m);
        for k = 1:m
            i=Ind_i(k); j=Ind_j(k); 
            Rij_est(:,:,k)=R_est(:,:,i)*(R_est(:,:,j)');
        end
        
        MSE = zeros(d,d,m);
        for j = 1:d
          MSE = MSE + bsxfun(@times,Rij_orig(:,j,:),Rij_est(:,j,:));
        end
        
        
        MSE_trace = zeros(1,m);
        
        for j = 1:d
            MSE_trace = MSE_trace + reshape(MSE(j,j,:),[1,m]);
        end
        
        
        mean_err = sqrt(mean(abs(1-MSE_trace/d)*0.5));

end