
    function lf = logfact(m)
        % compute the log factortial of m !
        if (m < 0)
            error('factorial should not be negative')
        elseif (m==0)
            lf = 0.;
            return
        elseif (m==1)
            lf = 0.;
            return
        end
        for k = 1:(m-1)
            if (k==1)
                lfm = 0;
            end
            lfmp1 = lfm + log(k+1);
            lfm = lfmp1;
                        
        end
        lf =lfmp1;
        
    end