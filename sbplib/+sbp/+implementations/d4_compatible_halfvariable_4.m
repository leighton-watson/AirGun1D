function [H, HI, D2, D4, e_1, e_m, M4, S2_1, S2_m, S3_1,...
    S3_m, S_1, S_m] = d4_compatible_halfvariable_4(m,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 4:de ordn. SBP Finita differens         %%%
    %%%                                         %%%
    %%% H           (Normen)                    %%%
    %%% D1=H^(-1)Q  (approx f?rsta derivatan)   %%%
    %%% D2          (approx andra derivatan)    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %m=20; %problemstorlek
    %h=1/(m-1);
    %h=1;
    
    BP = 6;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

%     c=ones(m,1);


    H=speye(m,m);
    H(1:4,1:4)=diag([17/48 59/48 43/48 49/48]);
    H(m-3:m,m-3:m)=rot90(diag([17/48 59/48 43/48 49/48]),2);
    H=H*h;
    HI=inv(H);
    HI = sparse(HI);



%     Q=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)-8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));
      e=ones(m,1);
%       Q=spdiags([e -8*e 0*e 8*e -e], -2:2, m, m)/12;
%     Q_U = [0 0.59e2 / 0.96e2 -0.1e1 / 0.12e2 -0.1e1 / 0.32e2; -0.59e2 / 0.96e2 0 0.59e2 / 0.96e2 0; 0.1e1 / 0.12e2 -0.59e2 / 0.96e2 0 0.59e2 / 0.96e2; 0.1e1 / 0.32e2 0 -0.59e2 / 0.96e2 0;];
%     Q(1:4,1:4)=Q_U;
%     Q(m-3:m,m-3:m)=rot90( -Q_U(1:4,1:4) ,2 );

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

%     D1=HI*(Q-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    M_U=[0.9e1 / 0.8e1 -0.59e2 / 0.48e2 0.1e1 / 0.12e2 0.1e1 / 0.48e2; -0.59e2 / 0.48e2 0.59e2 / 0.24e2 -0.59e2 / 0.48e2 0; 0.1e1 / 0.12e2 -0.59e2 / 0.48e2 0.55e2 / 0.24e2 -0.59e2 / 0.48e2; 0.1e1 / 0.48e2 0 -0.59e2 / 0.48e2 0.59e2 / 0.24e2;];
%     M=-(-1/12*diag(ones(m-2,1),2)+16/12*diag(ones(m-1,1),1)+16/12*diag(ones(m-1,1),-1)-1/12*diag(ones(m-2,1),-2)-30/12*diag(ones(m,1),0));
    M=-spdiags([-e 16*e -30*e 16*e -e], -2:2, m, m)/12;

    M(1:4,1:4)=M_U;

    M(m-3:m,m-3:m)=rot90(  M_U ,2 );
    M=M/h;

    S_U=[-0.11e2 / 0.6e1 3 -0.3e1 / 0.2e1 0.1e1 / 0.3e1;]/h;
    S_1=sparse(1,m);
    S_1(1:4)=S_U;
    S_m=sparse(1,m);
    S_m(m-3:m)=fliplr(-S_U);
    S_1 = S_1';
    S_m = S_m';


    M=sparse(m,m);
    e_1 = sparse(e_1);
    e_m = sparse(e_m);
    S_1 = sparse(S_1);
    S_m = sparse(S_m);


    scheme_width = 5;
    scheme_radius = (scheme_width-1)/2;
    r = (1+scheme_radius):(m-scheme_radius);

    function D2 = D2_fun(c)

        % ALTERNATIVES %%%%%%%%%%%%%
        % for i=4:m-3
        %     M(i,i-2:i+2)=[-c(i-1) / 0.6e1 + c(i-2) / 0.8e1 + c(i) / 0.8e1 -c(i-2) / 0.6e1 - c(i+1) / 0.6e1 - c(i-1) / 0.2e1 - c(i) / 0.2e1 c(i-2) / 0.24e2 + 0.5e1 / 0.6e1 * c(i-1) + 0.5e1 / 0.6e1 * c(i+1) + c(i+2) / 0.24e2 + 0.3e1 / 0.4e1 * c(i) -c(i-1) / 0.6e1 - c(i+2) / 0.6e1 - c(i) / 0.2e1 - c(i+1) / 0.2e1 -c(i+1) / 0.6e1 + c(i) / 0.8e1 + c(i+2) / 0.8e1;];
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for i=4:m-3
        %     M(i,i-2:i+2)= [
        %      c(i-2) / 0.8e1 - c(i-1) / 0.6e1 + c(i)   / 0.8e1                                   ,
        %     -c(i-2) / 0.6e1 - c(i-1) / 0.2e1 - c(i)   / 0.2e1 - c(i+1) / 0.6e1                  ,
        %      c(i-2) / 2.4e1 + c(i-1) / 1.2e0 + c(i) * 0.3/0.4 + c(i+1) / 1.2e0 + c(i+2) / 2.4e1 ,
        %                      -c(i-1) / 0.6e1 - c(i)   / 0.2e1 - c(i+1) / 0.2e1 - c(i+2) / 0.6e1 ,
        %                                        c(i)   / 0.8e1 - c(i+1) / 0.6e1 + c(i+2) / 0.8e1 ,
        %     ];
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mm2 =  c(r-2) / 0.8e1 - c(r-1) / 0.6e1 + c(r)   / 0.8e1                                  ;
        Mm1 = -c(r-2) / 0.6e1 - c(r-1) / 0.2e1 - c(r)   / 0.2e1 - c(r+1) / 0.6e1                 ;
        M0  =  c(r-2) / 2.4e1 + c(r-1) / 1.2e0 + c(r) * 0.3/0.4 + c(r+1) / 1.2e0 + c(r+2) / 2.4e1;
        Mp1 =                  -c(r-1) / 0.6e1 - c(r)   / 0.2e1 - c(r+1) / 0.2e1 - c(r+2) / 0.6e1;
        Mp2 =                                    c(r)   / 0.8e1 - c(r+1) / 0.6e1 + c(r+2) / 0.8e1;
        % printSize(Mm2);
        % scheme_radius
        % m
        M(r,:) = spdiags([Mm2 Mm1 M0 Mp1 Mp2],0:2*scheme_radius,length(r),m);
        % M(r,:) = spdiags([Mm2 Mm1 M0 Mp1 Mp2],(-2:2)+scheme_radius,M(r,:)); % This is slower
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %% Somthing is wrong here!!
        % Mm2 =  c(r-2) / 0.8e1 - c(r-1) / 0.6e1 + c(r)   / 0.8e1                                  ;
        % Mm1 = -c(r-2) / 0.6e1 - c(r-1) / 0.2e1 - c(r)   / 0.2e1 - c(r+1) / 0.6e1                 ;
        % M0  =  c(r-2) / 2.4e1 + c(r-1) / 1.2e0 + c(r) * 0.3/0.4 + c(r+1) / 1.2e0 + c(r+2) / 2.4e1;
        % Mp1 =                  -c(r-1) / 0.6e1 - c(r)   / 0.2e1 - c(r+1) / 0.2e1 - c(r+2) / 0.6e1;
        % Mp2 =                                    c(r)   / 0.8e1 - c(r+1) / 0.6e1 + c(r+2) / 0.8e1;
        % % printSize(M_diag_ind);
        % % Mdiags = [Mm2 Mm1 M0  Mp1 Mp2];
        % % printSize(Mdiags);
        % M(M_diag_ind) = [Mm2 Mm1 M0  Mp1 Mp2]; % This is slightly faster
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Kan man skriva det som en multiplikation av en 3-dim matris?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        M(1:6,1:6)=[0.12e2 / 0.17e2 * c(1) + 0.59e2 / 0.192e3 * c(2) + 0.27010400129e11 / 0.345067064608e12 * c(3) + 0.69462376031e11 / 0.2070402387648e13 * c(4) -0.59e2 / 0.68e2 * c(1) - 0.6025413881e10 / 0.21126554976e11 * c(3) - 0.537416663e9 / 0.7042184992e10 * c(4) 0.2e1 / 0.17e2 * c(1) - 0.59e2 / 0.192e3 * c(2) + 0.213318005e9 / 0.16049630912e11 * c(4) + 0.2083938599e10 / 0.8024815456e10 * c(3) 0.3e1 / 0.68e2 * c(1) - 0.1244724001e10 / 0.21126554976e11 * c(3) + 0.752806667e9 / 0.21126554976e11 * c(4) 0.49579087e8 / 0.10149031312e11 * c(3) - 0.49579087e8 / 0.10149031312e11 * c(4) -c(4) / 0.784e3 + c(3) / 0.784e3; -0.59e2 / 0.68e2 * c(1) - 0.6025413881e10 / 0.21126554976e11 * c(3) - 0.537416663e9 / 0.7042184992e10 * c(4) 0.3481e4 / 0.3264e4 * c(1) + 0.9258282831623875e16 / 0.7669235228057664e16 * c(3) + 0.236024329996203e15 / 0.1278205871342944e16 * c(4) -0.59e2 / 0.408e3 * c(1) - 0.29294615794607e14 / 0.29725717938208e14 * c(3) - 0.2944673881023e13 / 0.29725717938208e14 * c(4) -0.59e2 / 0.1088e4 * c(1) + 0.260297319232891e15 / 0.2556411742685888e16 * c(3) - 0.60834186813841e14 / 0.1278205871342944e16 * c(4) -0.1328188692663e13 / 0.37594290333616e14 * c(3) + 0.1328188692663e13 / 0.37594290333616e14 * c(4) -0.8673e4 / 0.2904112e7 * c(3) + 0.8673e4 / 0.2904112e7 * c(4); 0.2e1 / 0.17e2 * c(1) - 0.59e2 / 0.192e3 * c(2) + 0.213318005e9 / 0.16049630912e11 * c(4) + 0.2083938599e10 / 0.8024815456e10 * c(3) -0.59e2 / 0.408e3 * c(1) - 0.29294615794607e14 / 0.29725717938208e14 * c(3) - 0.2944673881023e13 / 0.29725717938208e14 * c(4) c(1) / 0.51e2 + 0.59e2 / 0.192e3 * c(2) + 0.13777050223300597e17 / 0.26218083221499456e17 * c(4) + 0.564461e6 / 0.13384296e8 * c(5) + 0.378288882302546512209e21 / 0.270764341349677687456e21 * c(3) c(1) / 0.136e3 - 0.125059e6 / 0.743572e6 * c(5) - 0.4836340090442187227e19 / 0.5525802884687299744e19 * c(3) - 0.17220493277981e14 / 0.89177153814624e14 * c(4) -0.10532412077335e14 / 0.42840005263888e14 * c(4) + 0.1613976761032884305e19 / 0.7963657098519931984e19 * c(3) + 0.564461e6 / 0.4461432e7 * c(5) -0.960119e6 / 0.1280713392e10 * c(4) - 0.3391e4 / 0.6692148e7 * c(5) + 0.33235054191e11 / 0.26452850508784e14 * c(3); 0.3e1 / 0.68e2 * c(1) - 0.1244724001e10 / 0.21126554976e11 * c(3) + 0.752806667e9 / 0.21126554976e11 * c(4) -0.59e2 / 0.1088e4 * c(1) + 0.260297319232891e15 / 0.2556411742685888e16 * c(3) - 0.60834186813841e14 / 0.1278205871342944e16 * c(4) c(1) / 0.136e3 - 0.125059e6 / 0.743572e6 * c(5) - 0.4836340090442187227e19 / 0.5525802884687299744e19 * c(3) - 0.17220493277981e14 / 0.89177153814624e14 * c(4) 0.3e1 / 0.1088e4 * c(1) + 0.507284006600757858213e21 / 0.475219048083107777984e21 * c(3) + 0.1869103e7 / 0.2230716e7 * c(5) + c(6) / 0.24e2 + 0.1950062198436997e16 / 0.3834617614028832e16 * c(4) -0.4959271814984644613e19 / 0.20965546238960637264e20 * c(3) - c(6) / 0.6e1 - 0.15998714909649e14 / 0.37594290333616e14 * c(4) - 0.375177e6 / 0.743572e6 * c(5) -0.368395e6 / 0.2230716e7 * c(5) + 0.752806667e9 / 0.539854092016e12 * c(3) + 0.1063649e7 / 0.8712336e7 * c(4) + c(6) / 0.8e1; 0.49579087e8 / 0.10149031312e11 * c(3) - 0.49579087e8 / 0.10149031312e11 * c(4) -0.1328188692663e13 / 0.37594290333616e14 * c(3) + 0.1328188692663e13 / 0.37594290333616e14 * c(4) -0.10532412077335e14 / 0.42840005263888e14 * c(4) + 0.1613976761032884305e19 / 0.7963657098519931984e19 * c(3) + 0.564461e6 / 0.4461432e7 * c(5) -0.4959271814984644613e19 / 0.20965546238960637264e20 * c(3) - c(6) / 0.6e1 - 0.15998714909649e14 / 0.37594290333616e14 * c(4) - 0.375177e6 / 0.743572e6 * c(5) 0.8386761355510099813e19 / 0.128413970713633903242e21 * c(3) + 0.2224717261773437e16 / 0.2763180339520776e16 * c(4) + 0.5e1 / 0.6e1 * c(6) + c(7) / 0.24e2 + 0.280535e6 / 0.371786e6 * c(5) -0.35039615e8 / 0.213452232e9 * c(4) - c(7) / 0.6e1 - 0.13091810925e11 / 0.13226425254392e14 * c(3) - 0.1118749e7 / 0.2230716e7 * c(5) - c(6) / 0.2e1; -c(4) / 0.784e3 + c(3) / 0.784e3 -0.8673e4 / 0.2904112e7 * c(3) + 0.8673e4 / 0.2904112e7 * c(4) -0.960119e6 / 0.1280713392e10 * c(4) - 0.3391e4 / 0.6692148e7 * c(5) + 0.33235054191e11 / 0.26452850508784e14 * c(3) -0.368395e6 / 0.2230716e7 * c(5) + 0.752806667e9 / 0.539854092016e12 * c(3) + 0.1063649e7 / 0.8712336e7 * c(4) + c(6) / 0.8e1 -0.35039615e8 / 0.213452232e9 * c(4) - c(7) / 0.6e1 - 0.13091810925e11 / 0.13226425254392e14 * c(3) - 0.1118749e7 / 0.2230716e7 * c(5) - c(6) / 0.2e1 0.3290636e7 / 0.80044587e8 * c(4) + 0.5580181e7 / 0.6692148e7 * c(5) + 0.5e1 / 0.6e1 * c(7) + c(8) / 0.24e2 + 0.660204843e9 / 0.13226425254392e14 * c(3) + 0.3e1 / 0.4e1 * c(6);];

        M(m-5:m,m-5:m)=[c(m-7) / 0.24e2 + 0.5e1 / 0.6e1 * c(m-6) + 0.5580181e7 / 0.6692148e7 * c(m-4) + 0.4887707739997e13 / 0.119037827289528e15 * c(m-3) + 0.3e1 / 0.4e1 * c(m-5) + 0.660204843e9 / 0.13226425254392e14 * c(m-2) + 0.660204843e9 / 0.13226425254392e14 * c(m-1) -c(m-6) / 0.6e1 - 0.1618585929605e13 / 0.9919818940794e13 * c(m-3) - c(m-5) / 0.2e1 - 0.1118749e7 / 0.2230716e7 * c(m-4) - 0.13091810925e11 / 0.13226425254392e14 * c(m-2) - 0.13091810925e11 / 0.13226425254392e14 * c(m-1) -0.368395e6 / 0.2230716e7 * c(m-4) + c(m-5) / 0.8e1 + 0.48866620889e11 / 0.404890569012e12 * c(m-3) + 0.752806667e9 / 0.539854092016e12 * c(m-2) + 0.752806667e9 / 0.539854092016e12 * c(m-1) -0.3391e4 / 0.6692148e7 * c(m-4) - 0.238797444493e12 / 0.119037827289528e15 * c(m-3) + 0.33235054191e11 / 0.26452850508784e14 * c(m-2) + 0.33235054191e11 / 0.26452850508784e14 * c(m-1) -0.8673e4 / 0.2904112e7 * c(m-2) - 0.8673e4 / 0.2904112e7 * c(m-1) + 0.8673e4 / 0.1452056e7 * c(m-3) -c(m-3) / 0.392e3 + c(m-2) / 0.784e3 + c(m-1) / 0.784e3; -c(m-6) / 0.6e1 - 0.1618585929605e13 / 0.9919818940794e13 * c(m-3) - c(m-5) / 0.2e1 - 0.1118749e7 / 0.2230716e7 * c(m-4) - 0.13091810925e11 / 0.13226425254392e14 * c(m-2) - 0.13091810925e11 / 0.13226425254392e14 * c(m-1) c(m-6) / 0.24e2 + 0.5e1 / 0.6e1 * c(m-5) + 0.3896014498639e13 / 0.4959909470397e13 * c(m-3) + 0.8386761355510099813e19 / 0.128413970713633903242e21 * c(m-2) + 0.280535e6 / 0.371786e6 * c(m-4) + 0.3360696339136261875e19 / 0.171218627618178537656e21 * c(m-1) -c(m-5) / 0.6e1 - 0.4959271814984644613e19 / 0.20965546238960637264e20 * c(m-2) - 0.375177e6 / 0.743572e6 * c(m-4) - 0.13425842714e11 / 0.33740880751e11 * c(m-3) - 0.193247108773400725e18 / 0.6988515412986879088e19 * c(m-1) -0.365281640980e12 / 0.1653303156799e13 * c(m-3) + 0.564461e6 / 0.4461432e7 * c(m-4) + 0.1613976761032884305e19 / 0.7963657098519931984e19 * c(m-2) - 0.198407225513315475e18 / 0.7963657098519931984e19 * c(m-1) -0.1328188692663e13 / 0.37594290333616e14 * c(m-2) + 0.2226377963775e13 / 0.37594290333616e14 * c(m-1) - 0.8673e4 / 0.363014e6 * c(m-3) c(m-3) / 0.49e2 + 0.49579087e8 / 0.10149031312e11 * c(m-2) - 0.256702175e9 / 0.10149031312e11 * c(m-1); -0.368395e6 / 0.2230716e7 * c(m-4) + c(m-5) / 0.8e1 + 0.48866620889e11 / 0.404890569012e12 * c(m-3) + 0.752806667e9 / 0.539854092016e12 * c(m-2) + 0.752806667e9 / 0.539854092016e12 * c(m-1) -c(m-5) / 0.6e1 - 0.4959271814984644613e19 / 0.20965546238960637264e20 * c(m-2) - 0.375177e6 / 0.743572e6 * c(m-4) - 0.13425842714e11 / 0.33740880751e11 * c(m-3) - 0.193247108773400725e18 / 0.6988515412986879088e19 * c(m-1) c(m-5) / 0.24e2 + 0.1869103e7 / 0.2230716e7 * c(m-4) + 0.507284006600757858213e21 / 0.475219048083107777984e21 * c(m-2) + 0.3e1 / 0.1088e4 * c(m) + 0.31688435395e11 / 0.67481761502e11 * c(m-3) + 0.27769176016102795561e20 / 0.712828572124661666976e21 * c(m-1) -0.125059e6 / 0.743572e6 * c(m-4) + c(m) / 0.136e3 - 0.23099342648e11 / 0.101222642253e12 * c(m-3) - 0.4836340090442187227e19 / 0.5525802884687299744e19 * c(m-2) + 0.193950157930938693e18 / 0.5525802884687299744e19 * c(m-1) 0.260297319232891e15 / 0.2556411742685888e16 * c(m-2) - 0.59e2 / 0.1088e4 * c(m) - 0.106641839640553e15 / 0.1278205871342944e16 * c(m-1) + 0.26019e5 / 0.726028e6 * c(m-3) -0.1244724001e10 / 0.21126554976e11 * c(m-2) + 0.3e1 / 0.68e2 * c(m) + 0.752806667e9 / 0.21126554976e11 * c(m-1); -0.3391e4 / 0.6692148e7 * c(m-4) - 0.238797444493e12 / 0.119037827289528e15 * c(m-3) + 0.33235054191e11 / 0.26452850508784e14 * c(m-2) + 0.33235054191e11 / 0.26452850508784e14 * c(m-1) -0.365281640980e12 / 0.1653303156799e13 * c(m-3) + 0.564461e6 / 0.4461432e7 * c(m-4) + 0.1613976761032884305e19 / 0.7963657098519931984e19 * c(m-2) - 0.198407225513315475e18 / 0.7963657098519931984e19 * c(m-1) -0.125059e6 / 0.743572e6 * c(m-4) + c(m) / 0.136e3 - 0.23099342648e11 / 0.101222642253e12 * c(m-3) - 0.4836340090442187227e19 / 0.5525802884687299744e19 * c(m-2) + 0.193950157930938693e18 / 0.5525802884687299744e19 * c(m-1) 0.564461e6 / 0.13384296e8 * c(m-4) + 0.470299699916357e15 / 0.952302618316224e15 * c(m-3) + 0.550597048646198778781e21 / 0.1624586048098066124736e22 * c(m-1) + c(m) / 0.51e2 + 0.378288882302546512209e21 / 0.270764341349677687456e21 * c(m-2) -0.59e2 / 0.408e3 * c(m) - 0.29294615794607e14 / 0.29725717938208e14 * c(m-2) - 0.2234477713167e13 / 0.29725717938208e14 * c(m-1) - 0.8673e4 / 0.363014e6 * c(m-3) -0.59e2 / 0.3136e4 * c(m-3) - 0.13249937023e11 / 0.48148892736e11 * c(m-1) + 0.2e1 / 0.17e2 * c(m) + 0.2083938599e10 / 0.8024815456e10 * c(m-2); -0.8673e4 / 0.2904112e7 * c(m-2) - 0.8673e4 / 0.2904112e7 * c(m-1) + 0.8673e4 / 0.1452056e7 * c(m-3) -0.1328188692663e13 / 0.37594290333616e14 * c(m-2) + 0.2226377963775e13 / 0.37594290333616e14 * c(m-1) - 0.8673e4 / 0.363014e6 * c(m-3) 0.260297319232891e15 / 0.2556411742685888e16 * c(m-2) - 0.59e2 / 0.1088e4 * c(m) - 0.106641839640553e15 / 0.1278205871342944e16 * c(m-1) + 0.26019e5 / 0.726028e6 * c(m-3) -0.59e2 / 0.408e3 * c(m) - 0.29294615794607e14 / 0.29725717938208e14 * c(m-2) - 0.2234477713167e13 / 0.29725717938208e14 * c(m-1) - 0.8673e4 / 0.363014e6 * c(m-3) 0.9258282831623875e16 / 0.7669235228057664e16 * c(m-2) + 0.3481e4 / 0.3264e4 * c(m) + 0.228389721191751e15 / 0.1278205871342944e16 * c(m-1) + 0.8673e4 / 0.1452056e7 * c(m-3) -0.6025413881e10 / 0.21126554976e11 * c(m-2) - 0.59e2 / 0.68e2 * c(m) - 0.537416663e9 / 0.7042184992e10 * c(m-1); -c(m-3) / 0.392e3 + c(m-2) / 0.784e3 + c(m-1) / 0.784e3 c(m-3) / 0.49e2 + 0.49579087e8 / 0.10149031312e11 * c(m-2) - 0.256702175e9 / 0.10149031312e11 * c(m-1) -0.1244724001e10 / 0.21126554976e11 * c(m-2) + 0.3e1 / 0.68e2 * c(m) + 0.752806667e9 / 0.21126554976e11 * c(m-1) -0.59e2 / 0.3136e4 * c(m-3) - 0.13249937023e11 / 0.48148892736e11 * c(m-1) + 0.2e1 / 0.17e2 * c(m) + 0.2083938599e10 / 0.8024815456e10 * c(m-2) -0.6025413881e10 / 0.21126554976e11 * c(m-2) - 0.59e2 / 0.68e2 * c(m) - 0.537416663e9 / 0.7042184992e10 * c(m-1) 0.3e1 / 0.3136e4 * c(m-3) + 0.27010400129e11 / 0.345067064608e12 * c(m-2) + 0.234566387291e12 / 0.690134129216e12 * c(m-1) + 0.12e2 / 0.17e2 * c(m);];

        M=M/h;
        D2=HI*(-M-c(1)*e_1*S_1'+c(m)*e_m*S_m');
    end
    D2 = @D2_fun;


    S2_U=[2 -5 4 -1;]/h^2;
    S2_1=sparse(1,m);
    S2_1(1:4)=S2_U;
    S2_m=sparse(1,m);
    S2_m(m-3:m)=fliplr(S2_U);
    S2_1 = S2_1';
    S2_m = S2_m';

    m3=-1/6;m2=2;m1=-13/2;m0=28/3;
%     M4=m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);
    stencil = [m3,m2,m1,m0,m1,m2,m3];
    d = (length(stencil)-1)/2;
    diags = -d:d;
    M4 = stripeMatrix(stencil, diags, m);

    %M4=(-1/6*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3) ) + 2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2)) -13/2*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)) + 28/3*diag(ones(m,1),0));

    M4_U=[0.5762947e7 / 0.2316384e7 -0.6374287e7 / 0.1158192e7 0.573947e6 / 0.165456e6 -0.124637e6 / 0.289548e6 0.67979e5 / 0.2316384e7 -0.60257e5 / 0.1158192e7; -0.6374287e7 / 0.1158192e7 0.30392389e8 / 0.2316384e7 -0.2735053e7 / 0.289548e6 0.273109e6 / 0.165456e6 0.83767e5 / 0.1158192e7 0.245549e6 / 0.2316384e7; 0.573947e6 / 0.165456e6 -0.2735053e7 / 0.289548e6 0.5266855e7 / 0.579096e6 -0.1099715e7 / 0.289548e6 0.869293e6 / 0.1158192e7 -0.10195e5 / 0.144774e6; -0.124637e6 / 0.289548e6 0.273109e6 / 0.165456e6 -0.1099715e7 / 0.289548e6 0.3259225e7 / 0.579096e6 -0.324229e6 / 0.72387e5 0.1847891e7 / 0.1158192e7; 0.67979e5 / 0.2316384e7 0.83767e5 / 0.1158192e7 0.869293e6 / 0.1158192e7 -0.324229e6 / 0.72387e5 0.2626501e7 / 0.330912e6 -0.7115491e7 / 0.1158192e7; -0.60257e5 / 0.1158192e7 0.245549e6 / 0.2316384e7 -0.10195e5 / 0.144774e6 0.1847891e7 / 0.1158192e7 -0.7115491e7 / 0.1158192e7 0.21383077e8 / 0.2316384e7;];

    M4(1:6,1:6)=M4_U;

    M4(m-5:m,m-5:m)=rot90(  M4_U ,2 );
    M4=M4/h^3;

    S3_U=[-1 3 -3 1;]/h^3;
    S3_1=sparse(1,m);
    S3_1(1:4)=S3_U;
    S3_m=sparse(1,m);
    S3_m(m-3:m)=fliplr(-S3_U);
    S3_1 = S3_1';
    S3_m = S3_m';

    D4=HI*(M4-e_1*S3_1'+e_m*S3_m'  + S_1*S2_1'-S_m*S2_m');





end