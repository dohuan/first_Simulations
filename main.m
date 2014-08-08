% =========================================================================
%       ---------------------------------------------------------
%        Test the effect of Tetha and its convergance properties
%       ---------------------------------------------------------
% Authour : Mahdi Jadaliha
% e-mail  : jadaliha@gmail.com
% =========================================================================
%tic  
close all
clear all
clc
 
%----- initialize system configuration-------------------------------------
global globaloption  sys_parameter
[ globaloption , sys_parameter ] = Configuration();
%----- generate or load the field------------------------------------------
option = globaloption;
[true_field, VehicleStates, AgumentedData, ~]  =  ...                       % [field, VehicleStates, AgumentedData, numberofpossibilities] = Generate_new_field(savetofile)
    Generate_new_field();                                                   % the data will be saved in 'savetofile'

%--------------------------------------------------------------------------
nx = size(option.X_mesh,2) ;                                                 % number of X grid
ny = size(option.Y_mesh,2) ;                                                % number of Y grid
nt = size(option.T_mesh,2) ;                                                % number of time steps
ntheta= size(option.hyperparameters_possibilities,1);                       % number of possibilities for $theta$
n = size(option.grids,1);                                                   % number of spatial sites
N = option.agentnumbers;
x_width = option.X_mesh(end) - option.X_mesh(1) + option.finegridsize;
y_width = option.Y_mesh(end) - option.Y_mesh(1) + option.finegridsize;

%----- Construct precision matrix Qz---------------------------------------
muz_theta = zeros(n,ntheta);
Sigmaz_theta = zeros(n^2,ntheta);
for indtheta=1:ntheta 
    tmp1 = prec_mat_5by5(option.hyperparameters_possibilities(indtheta,1)...
        ,option.hyperparameters_possibilities(indtheta,2),nx,ny);
    Sigmaz_theta(:,indtheta) = reshape(tmp1^(-1),[],1);
end 


alpha = 0.6;

f_theta = log(option.hyperparameters_possibilities(:,end));
itmp = (1:2*N);
jtmp = kron((1:N)',[1;1]);
mux = zeros(2*N,1);
Sigmax = eye(2*N)*1000;
mux_ = mux;
Sigmax_ = Sigmax;

for t=option.T_mesh 
    t
    %tic
    xtilda = reshape(option.grids(...
        AgumentedData(1,t).possible_q.measuredposition,:)',[],1);
    ztilda = AgumentedData(t).y;

    
    qt= AgumentedData(1,t).possible_q.support_qt{1,1}';
    for indj = 2:N
        qtN = AgumentedData(1,t).possible_q.support_qt{1,indj}';
        qt= [kron(ones(1,size(qtN,2)),qt);...
            kron(qtN,ones(1,size(qt,2)))];
    end
    Sigma_e = eye(2*size(qt,1))* option.vehicle.ObservationNoise;


    
    numberofpossibleqt = size(qt,2);
    Sigma_epsilon = eye(size(qt,1))* option.s_e2;
    f_qtandtheta = zeros(numberofpossibleqt,ntheta);
    muz_qandtheta =    zeros(n, 1);
    Sigmaz_qandtheta = zeros(n^2,1);
    %tic
    for indq = 1:numberofpossibleqt
        q = qt(:,indq);
        H = sparse(1:N,q',ones(1,N),N,n);                                   % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
        x = reshape(option.grids(q,:)',[],1);
        f_qt  = logmvnpdf2(x,mux_,Sigmax_,option);
        f_qtildaGqt  = logmvnpdf2(xtilda,x,Sigma_e,option);
        for indtheta = 1:ntheta
            muztheta = muz_theta(:,indtheta);
            Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
            muztildatheta    = H * muztheta;                                % $mu_{tilde{zt}|theta,D{t-1},\q{t}}  =  H_{qt} mu_{z| theta,D{t-1}}$            
            Sigmaztildatheta = Sigma_epsilon + H * Sigmaztheta * H';
            f_ztilda = logmvnpdf(ztilda,muztildatheta,Sigmaztildatheta);
            % approximate the distribution of $\q{t},\theta |\D{0}{t}$           
            f_qtandtheta(indq,indtheta) = ...
                f_qtildaGqt  + f_theta(indtheta) + f_qt + f_ztilda;         % pi(q{t},theta|D{t})
            

        end
    end
    %toc
    tmp_c = log(sum(sum(exp(f_qtandtheta))));                               % Here we find normalization factor $tmp_c$
%     if (abs(tmp_c)>100)
%         pause
%     end
    f_qtandtheta = f_qtandtheta - tmp_c;                                    % After this line $f_qtandtheta$ is normalized, 
    pi_qtandtheta = exp(f_qtandtheta);                                      % $sum(pi_qtandtheta) = 1$.
    
    pi_theta = sum(pi_qtandtheta,1);                
    pi_theta = pi_theta/sum(pi_theta);                                      % Compensate computational error
    f_theta = log(pi_theta);
    pi_q = sum(pi_qtandtheta,2);
    
    
    
    
    % For the toures correction on the sampling positions measured  close
    % to the borders
    mux = zeros(2*N,1);
    xxx = zeros(2*N,numberofpossibleqt);
    for indj = 1:N
        xx = option.grids(qt(indj,:),:);
        median_xx = median(xx);
        tmp1 = find((xx(:,1) - median_xx(1))> x_width/2);
        xx(tmp1,1) = xx(tmp1,1)- x_width;
        tmp2 = find((xx(:,1) - median_xx(1))<-x_width/2);
        xx(tmp2,1) = xx(tmp2,1)+ x_width;
        tmp1 = find((xx(:,2) - median_xx(2))> y_width/2);
        xx(tmp1,2) = xx(tmp1,2)- y_width;
        tmp2 = find((xx(:,2) - median_xx(2))<-y_width/2);
        xx(tmp2,2) = xx(tmp2,2)+ y_width;
        
        xxx((indj*2)-1:indj*2,:) = xx';
    end
    mux = (xxx * pi_q);
    Sigmax = reshape(...
        VectorSquare(xxx - mux * ones(1,numberofpossibleqt)) * pi_q,...
        2*N,2*N);



    
    
    
    
    
    
    
    tic
    for indtheta = 1:ntheta
        pi_qtGtheta = (pi_qtandtheta(:,indtheta)/pi_theta(indtheta));
        muz_theta_temp = zeros(n,1);
        Sigmaz_theta_temp = zeros(n^2,1);
        for indq = 1:numberofpossibleqt
            q = qt(:,indq);
            H = sparse(1:N,q',ones(1,N),N,n);                               % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
            muztheta = muz_theta(:,indtheta);
            Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
            Sigmaztheta = 0.5 * (Sigmaztheta+Sigmaztheta');
            muztildatheta    = H * muztheta;
            tmp1 = Sigmaztheta * H';
            Sigmaztildatheta = Sigma_epsilon + H * tmp1;
            InvSigmaztildatheta = Sigmaztildatheta^(-1);
 
            muz_qandtheta    = ...
                muztheta + tmp1 * ...
                InvSigmaztildatheta * (ztilda - muztildatheta);
            Sigmaz_qandtheta = ...
                reshape(Sigmaztheta - tmp1 * ...
                InvSigmaztildatheta * tmp1',[],1);
            
            muz_theta_temp = muz_theta_temp + ...
                muz_qandtheta * pi_qtGtheta(indq);

            Sigmaz_theta_temp = Sigmaz_theta_temp + ...
                (Sigmaz_qandtheta + VectorSquare(muz_qandtheta))*...
                pi_qtGtheta(indq) ;

        end
        muz_theta(:,indtheta)    = muz_theta_temp;
        Sigmaz_theta(:,indtheta) = ...
            Sigmaz_theta_temp - VectorSquare(muz_theta_temp);
    end
    toc
    posterior(t).muz = muz_theta * pi_theta';
    posterior(t).varz= ...
        (Sigmaz_theta(1:n+1:end,:)+muz_theta.*muz_theta)* pi_theta' ...
        - posterior(t).muz.* posterior(t).muz;
    posterior(t).mux = mux;
    posterior(t).Sigmax = Sigmax;

    % predict next sampling position
    u = VehicleStates(t).u;
    phi = VehicleStates(1,t).h';%
    stmp = reshape([cos(phi),sin(phi)]',[],1);
    F = sparse(itmp,jtmp,stmp,2*N,N);
    %mux_    = mux + F * u;
    mux_    = alpha*mux;
    Sigma_w = eye(2*size(qt,1))* option.vehicle.ModelUncertanity;
    Sigmax_ = Sigmax + Sigma_w;  
    
    toc
end
% 
toc
% 
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %-----------------------    True Position   -------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% nx = size(option.X_mesh,2) ;                                                % number of X grid
% ny = size(option.Y_mesh,2) ;                                                % number of Y grid
% nt = size(option.T_mesh,2) ;                                                % number of time steps
% ntheta= size(option.hyperparameters_possibilities,1);                       % number of possibilities for $theta$
% n = size(option.grids,1);                                                   % number of spatial sites
% N = option.agentnumbers;
% x_width = option.X_mesh(end) - option.X_mesh(1) + option.finegridsize;
% y_width = option.Y_mesh(end) - option.Y_mesh(1) + option.finegridsize;
% 
% %----- Construct precision matrix Qz---------------------------------------
% muz_theta = zeros(n,ntheta);
% Sigmaz_theta = zeros(n^2,ntheta);
% for indtheta=1:ntheta 
%     tmp1 = prec_mat_5by5(option.hyperparameters_possibilities(indtheta,1)...
%         ,option.hyperparameters_possibilities(indtheta,2),nx,ny);
%     Sigmaz_theta(:,indtheta) = reshape(tmp1^(-1),[],1);
% end 
% 
% 
% 
% f_theta = log(option.hyperparameters_possibilities(:,end));
% itmp = (1:2*N);
% jtmp = kron((1:N)',[1;1]);
% mux = zeros(2*N,1);
% Sigmax = eye(2*N)*1000;
% mux_ = mux;
% Sigmax_ = Sigmax;
% 
% for t=option.T_mesh 
%     t
%     xtilda = reshape(option.grids(...
%         AgumentedData(1,t).possible_q.measuredposition,:)',[],1);
%     ztilda = AgumentedData(t).y;
% 
%     
%     qt= AgumentedData(1,t).possible_q.true_q';
%    
%     Sigma_e = eye(2*size(qt,1))* option.vehicle.ObservationNoise;
% 
% 
%     
%     numberofpossibleqt = size(qt,2);
%     Sigma_epsilon = eye(size(qt,1))* option.s_e2;
%     f_qtandtheta = zeros(numberofpossibleqt,ntheta);
%     muz_qandtheta =    zeros(n, 1);
%     Sigmaz_qandtheta = zeros(n^2,1);
%     for indq = 1:numberofpossibleqt
%         q = qt(:,indq);
%         H = sparse(1:N,q',ones(1,N),N,n);                                   % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
%         x = reshape(option.grids(q,:)',[],1);
%         f_qt  = logmvnpdf2(x,mux_,Sigmax_,option);
%         f_qtildaGqt  = logmvnpdf2(xtilda,x,Sigma_e,option);
%         for indtheta = 1:ntheta
%             muztheta = muz_theta(:,indtheta);
%             Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
%             muztildatheta    = H * muztheta;                                % $mu_{tilde{zt}|theta,D{t-1},\q{t}}  =  H_{qt} mu_{z| theta,D{t-1}}$            
%             Sigmaztildatheta = Sigma_epsilon + H * Sigmaztheta * H';
%             f_ztilda = logmvnpdf(ztilda,muztildatheta,Sigmaztildatheta);
%             % approximate the distribution of $\q{t},\theta |\D{0}{t}$           
%             f_qtandtheta(indq,indtheta) = ...
%                 f_qtildaGqt  + f_theta(indtheta) + f_qt + f_ztilda;                % pi(q{t},theta|D{t})
%             
% 
%         end
%     end
%     tmp_c = log(sum(sum(exp(f_qtandtheta))));
% %     if (abs(tmp_c)>100)
% %         pause
% %     end
%     f_qtandtheta = f_qtandtheta - tmp_c;
%     pi_qtandtheta = exp(f_qtandtheta);
%     
%     pi_theta = sum(pi_qtandtheta,1);
%     f_theta = log(pi_theta);
%     pi_q = sum(pi_qtandtheta,2);
%     
%     
%     
%     
%     % For the toures correction on the sampling positions measured  close
%     % to the borders
%     mux = zeros(2*N,1);
%     xxx = zeros(2*N,numberofpossibleqt);
%     for indj = 1:N
%         xx = option.grids(qt(indj,:),:);
%         
%         xxx((indj*2)-1:indj*2,:) = xx';
%     end
%     mux = (xxx * pi_q);
%     Sigmax = reshape(...
%         VectorSquare(xxx - mux * ones(1,numberofpossibleqt)) * pi_q,...
%         2*N,2*N);
% 
%     for indtheta = 1:ntheta
%         pi_qtGtheta = (pi_qtandtheta(:,indtheta)/pi_theta(indtheta));
%         muz_theta_temp = zeros(n,1);
%         Sigmaz_theta_temp = zeros(n^2,1);
%         for indq = 1:numberofpossibleqt
%             q = qt(:,indq);
%             H = sparse(1:N,q',ones(1,N),N,n);                               % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
%             muztheta = muz_theta(:,indtheta);
%             Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
%             Sigmaztheta = 0.5 * (Sigmaztheta+Sigmaztheta');
%             muztildatheta    = H * muztheta;
%             tmp1 = Sigmaztheta * H';
%             Sigmaztildatheta = Sigma_epsilon + H * tmp1;
%             InvSigmaztildatheta = Sigmaztildatheta^(-1);
%  
%             muz_qandtheta    = ...
%                 muztheta + tmp1 * ...
%                 InvSigmaztildatheta * (ztilda - muztildatheta);
%             Sigmaz_qandtheta = ...
%                 reshape(Sigmaztheta - tmp1 * ...
%                 InvSigmaztildatheta * tmp1',[],1);
%             
%             muz_theta_temp = muz_theta_temp + ...
%                 muz_qandtheta * pi_qtGtheta(indq);
% 
%             Sigmaz_theta_temp = Sigmaz_theta_temp + ...
%                 (Sigmaz_qandtheta + VectorSquare(muz_qandtheta))*...
%                 pi_qtGtheta(indq) ;
% 
%         end
%         muz_theta(:,indtheta)    = muz_theta_temp;
%         Sigmaz_theta(:,indtheta) = ...
%             Sigmaz_theta_temp - VectorSquare(muz_theta_temp);
%     end
%     
%     posteriorTrueQ(t).muz = muz_theta * pi_theta';
%     posteriorTrueQ(t).varz= ...
%         (Sigmaz_theta(1:n+1:end,:)+muz_theta.*muz_theta)* pi_theta' ...
%         - posterior(t).muz.* posterior(t).muz;
%     posteriorTrueQ(t).mux = mux;
%     posteriorTrueQ(t).Sigmax = Sigmax;
% 
%     % predict next sampling position
%     u = VehicleStates(t).u;
%     phi = VehicleStates(1,t).h';
%     stmp = reshape([cos(phi),sin(phi)]',[],1);
%     F = sparse(itmp,jtmp,stmp,2*N,N);
%     mux_    = mux + F * u;
%     Sigma_w = eye(2*size(qt,1))* option.vehicle.ModelUncertanity;
%     Sigmax_ = Sigmax + Sigma_w;  
%     
%     
% end
% 
% 
% 
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %-----------------------   Noisy Position   -------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% nx = size(option.X_mesh,2) ;                                                % number of X grid
% ny = size(option.Y_mesh,2) ;                                                % number of Y grid
% nt = size(option.T_mesh,2) ;                                                % number of time steps
% ntheta= size(option.hyperparameters_possibilities,1);                       % number of possibilities for $theta$
% n = size(option.grids,1);                                                   % number of spatial sites
% N = option.agentnumbers;
% x_width = option.X_mesh(end) - option.X_mesh(1) + option.finegridsize;
% y_width = option.Y_mesh(end) - option.Y_mesh(1) + option.finegridsize;
% 
% %----- Construct precision matrix Qz---------------------------------------
% muz_theta = zeros(n,ntheta);
% Sigmaz_theta = zeros(n^2,ntheta);
% for indtheta=1:ntheta 
%     tmp1 = prec_mat_5by5(option.hyperparameters_possibilities(indtheta,1)...
%         ,option.hyperparameters_possibilities(indtheta,2),nx,ny);
%     Sigmaz_theta(:,indtheta) = reshape(tmp1^(-1),[],1);
% end 
% 
% 
% 
% f_theta = log(option.hyperparameters_possibilities(:,end));
% itmp = (1:2*N);
% jtmp = kron((1:N)',[1;1]);
% mux = zeros(2*N,1);
% Sigmax = eye(2*N)*1000;
% mux_ = mux;
% Sigmax_ = Sigmax;
% 
% for t=option.T_mesh 
%     t
%     xtilda = reshape(option.grids(...
%         AgumentedData(1,t).possible_q.measuredposition,:)',[],1);
%     ztilda = AgumentedData(t).y;
% 
%     
%     qt= AgumentedData(1,t).possible_q.measuredposition';
%    
%     Sigma_e = eye(2*size(qt,1))* option.vehicle.ObservationNoise;
% 
% 
%     
%     numberofpossibleqt = size(qt,2);
%     Sigma_epsilon = eye(size(qt,1))* option.s_e2;
%     f_qtandtheta = zeros(numberofpossibleqt,ntheta);
%     muz_qandtheta =    zeros(n, 1);
%     Sigmaz_qandtheta = zeros(n^2,1);
%     for indq = 1:numberofpossibleqt
%         q = qt(:,indq);
%         H = sparse(1:N,q',ones(1,N),N,n);                                   % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
%         x = reshape(option.grids(q,:)',[],1);
%         f_qt  = logmvnpdf2(x,mux_,Sigmax_,option);
%         f_qtildaGqt  = logmvnpdf2(xtilda,x,Sigma_e,option);
%         for indtheta = 1:ntheta
%             muztheta = muz_theta(:,indtheta);
%             Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
%             muztildatheta    = H * muztheta;                                % $mu_{tilde{zt}|theta,D{t-1},\q{t}}  =  H_{qt} mu_{z| theta,D{t-1}}$            
%             Sigmaztildatheta = Sigma_epsilon + H * Sigmaztheta * H';
%             f_ztilda = logmvnpdf(ztilda,muztildatheta,Sigmaztildatheta);
%             % approximate the distribution of $\q{t},\theta |\D{0}{t}$           
%             f_qtandtheta(indq,indtheta) = ...
%                 f_qtildaGqt  + f_theta(indtheta) + f_qt + f_ztilda;                % pi(q{t},theta|D{t})
%             
% 
%         end
%     end
%     tmp_c = log(sum(sum(exp(f_qtandtheta))));
% %     if (abs(tmp_c)>100)
% %         pause
% %     end
%     f_qtandtheta = f_qtandtheta - tmp_c;
%     pi_qtandtheta = exp(f_qtandtheta);
%     
%     pi_theta = sum(pi_qtandtheta,1);
%     f_theta = log(pi_theta);
%     pi_q = sum(pi_qtandtheta,2);
%     
%     
%     
%     
%     % For the toures correction on the sampling positions measured  close
%     % to the borders
%     mux = zeros(2*N,1);
%     xxx = zeros(2*N,numberofpossibleqt);
%     for indj = 1:N
%         xx = option.grids(qt(indj,:),:);
% 
%         xxx((indj*2)-1:indj*2,:) = xx';
%     end
%     mux = (xxx * pi_q);
%     Sigmax = reshape(...
%         VectorSquare(xxx - mux * ones(1,numberofpossibleqt)) * pi_q,...
%         2*N,2*N);
% 
%     for indtheta = 1:ntheta
%         pi_qtGtheta = (pi_qtandtheta(:,indtheta)/pi_theta(indtheta));
%         muz_theta_temp = zeros(n,1);
%         Sigmaz_theta_temp = zeros(n^2,1);
%         for indq = 1:numberofpossibleqt
%             q = qt(:,indq);
%             H = sparse(1:N,q',ones(1,N),N,n);                               % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
%             muztheta = muz_theta(:,indtheta);
%             Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
%             Sigmaztheta = 0.5 * (Sigmaztheta+Sigmaztheta');
%             muztildatheta    = H * muztheta;
%             tmp1 = Sigmaztheta * H';
%             Sigmaztildatheta = Sigma_epsilon + H * tmp1;
%             InvSigmaztildatheta = Sigmaztildatheta^(-1);
%  
%             muz_qandtheta    = ...
%                 muztheta + tmp1 * ...
%                 InvSigmaztildatheta * (ztilda - muztildatheta);
%             Sigmaz_qandtheta = ...
%                 reshape(Sigmaztheta - tmp1 * ...
%                 InvSigmaztildatheta * tmp1',[],1);
%             
%             muz_theta_temp = muz_theta_temp + ...
%                 muz_qandtheta * pi_qtGtheta(indq);
% 
%             Sigmaz_theta_temp = Sigmaz_theta_temp + ...
%                 (Sigmaz_qandtheta + VectorSquare(muz_qandtheta))*...
%                 pi_qtGtheta(indq) ;
% 
%         end
%         muz_theta(:,indtheta)    = muz_theta_temp;
%         Sigmaz_theta(:,indtheta) = ...
%             Sigmaz_theta_temp - VectorSquare(muz_theta_temp);
%     end
%     
%     posteriorNoisyQ(t).muz = muz_theta * pi_theta';
%     posteriorNoisyQ(t).varz= ...
%         (Sigmaz_theta(1:n+1:end,:)+muz_theta.*muz_theta)* pi_theta' ...
%         - posterior(t).muz.* posterior(t).muz;
%     posteriorNoisyQ(t).mux = mux;
%     posteriorNoisyQ(t).Sigmax = Sigmax;
% 
%     % predict next sampling position
%     u = VehicleStates(t).u;
%     phi = VehicleStates(1,t).h';
%     stmp = reshape([cos(phi),sin(phi)]',[],1);
%     F = sparse(itmp,jtmp,stmp,2*N,N);
%     mux_    = mux + F * u;
%     Sigma_w = eye(2*size(qt,1))* option.vehicle.ModelUncertanity;
%     Sigmax_ = Sigmax + Sigma_w;  
%     
%     
% end

save('test13alien.mat', ...
    'posterior', 'posteriorTrueQ', 'posteriorNoisyQ',...
    'AgumentedData', 'true_field', ...
    'option', 'sys_parameter')
