function [H, HI, D1, D2, e_1, e_m, M, Q, S_1, S_m] = d2_blocknorm_6(m,h)
    
    BP = 6;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    H_U=[0.8489084265971e13 / 0.45952647390720e14 0.24636450459943e14 / 0.98469958694400e14 -0.2796787072531e13 / 0.12308744836800e14 0.2793599068823e13 / 0.14360202309600e14 -0.66344569931569e14 / 0.689289710860800e15 0.3784697867191e13 / 0.137857942172160e15; 0.24636450459943e14 / 0.98469958694400e14 0.27815394775103e14 / 0.19693991738880e14 -0.445601472229e12 / 0.861612138576e12 0.3896159037731e13 / 0.17232242771520e14 -0.866505556741e12 / 0.27571588434432e14 -0.25625418493681e14 / 0.689289710860800e15; -0.2796787072531e13 / 0.12308744836800e14 -0.445601472229e12 / 0.861612138576e12 0.31409405327129e14 / 0.17232242771520e14 -0.1595539040819e13 / 0.3446448554304e13 0.2651608170899e13 / 0.17232242771520e14 0.1434714163381e13 / 0.43080606928800e14; 0.2793599068823e13 / 0.14360202309600e14 0.3896159037731e13 / 0.17232242771520e14 -0.1595539040819e13 / 0.3446448554304e13 0.6984350202787e13 / 0.5744080923840e13 -0.62662743973e11 / 0.861612138576e12 -0.435331581619e12 / 0.12308744836800e14; -0.66344569931569e14 / 0.689289710860800e15 -0.866505556741e12 / 0.27571588434432e14 0.2651608170899e13 / 0.17232242771520e14 -0.62662743973e11 / 0.861612138576e12 0.20320736807807e14 / 0.19693991738880e14 0.1368363924007e13 / 0.98469958694400e14; 0.3784697867191e13 / 0.137857942172160e15 -0.25625418493681e14 / 0.689289710860800e15 0.1434714163381e13 / 0.43080606928800e14 -0.435331581619e12 / 0.12308744836800e14 0.1368363924007e13 / 0.98469958694400e14 0.27414523542149e14 / 0.27571588434432e14;];


    H=speye(m);
    H(1:6,1:6)=H_U;
    H(m-5:m,m-5:m)=rot90( H_U(1:6,1:6) ,2 );
    H=H*h;
    HI=inv(H);

    diags   = -3:3;
    stencil = [-1/60,3/20,-3/4,0,3/4,-3/20,1/60];
    Q = stripeMatrix(stencil, diags, m);
    
%     Q=(1/60*diag(ones(m-3,1),3)-9/60*diag(ones(m-2,1),2)+45/60*diag(ones(m-1,1),1)-45/60*diag(ones(m-1,1),-1)+9/60*diag(ones(m-2,1),-2)-1/60*diag(ones(m-3,1),-3));

    Q_U = [-0.1e1 / 0.2e1 0.151864337282617e15 / 0.172322427715200e15 -0.251539972254817e15 / 0.344644855430400e15 0.61230525943549e14 / 0.114881618476800e15 -0.80987306509439e14 / 0.344644855430400e15 0.697178163343e12 / 0.13785794217216e14; -0.151864337282617e15 / 0.172322427715200e15 0 0.12350422095979e14 / 0.7658774565120e13 -0.78802251164141e14 / 0.68928971086080e14 0.4229407848431e13 / 0.7658774565120e13 -0.5372490790279e13 / 0.38293872825600e14; 0.251539972254817e15 / 0.344644855430400e15 -0.12350422095979e14 / 0.7658774565120e13 0 0.2217674201683e13 / 0.1723224277152e13 -0.13219134462287e14 / 0.22976323695360e14 0.19660399553981e14 / 0.114881618476800e15; -0.61230525943549e14 / 0.114881618476800e15 0.78802251164141e14 / 0.68928971086080e14 -0.2217674201683e13 / 0.1723224277152e13 0 0.62307836637379e14 / 0.68928971086080e14 -0.84068101764193e14 / 0.344644855430400e15; 0.80987306509439e14 / 0.344644855430400e15 -0.4229407848431e13 / 0.7658774565120e13 0.13219134462287e14 / 0.22976323695360e14 -0.62307836637379e14 / 0.68928971086080e14 0 0.44756810052211e14 / 0.57440809238400e14; -0.697178163343e12 / 0.13785794217216e14 0.5372490790279e13 / 0.38293872825600e14 -0.19660399553981e14 / 0.114881618476800e15 0.84068101764193e14 / 0.344644855430400e15 -0.44756810052211e14 / 0.57440809238400e14 0;];

    Q(1:6,1:6)=Q_U;
    Q(m-5:m,m-5:m)=rot90( -Q_U(1:6,1:6) ,2 );

    D1=H\Q;

    M_U=[0.960901171090739e15 / 0.689289710860800e15 -0.502032138770899e15 / 0.229763236953600e15 0.493085196645929e15 / 0.344644855430400e15 -0.329491854944251e15 / 0.344644855430400e15 0.89541920186441e14 / 0.229763236953600e15 -0.50617198740721e14 / 0.689289710860800e15; -0.100483015499831e15 / 0.45952647390720e14 0.807564929223191e15 / 0.137857942172160e15 -0.415779274818991e15 / 0.68928971086080e14 0.80693719872887e14 / 0.22976323695360e14 -0.196663473955997e15 / 0.137857942172160e15 0.37943821632959e14 / 0.137857942172160e15; 0.99938177941669e14 / 0.68928971086080e14 -0.84419552767043e14 / 0.13785794217216e14 0.106922123424097e15 / 0.11488161847680e14 -0.223356054245897e15 / 0.34464485543040e14 0.157526160982357e15 / 0.68928971086080e14 -0.10062402380533e14 / 0.22976323695360e14; -0.68310884976863e14 / 0.68928971086080e14 0.17038649985979e14 / 0.4595264739072e13 -0.231397767539273e15 / 0.34464485543040e14 0.232669188399619e15 / 0.34464485543040e14 -0.1657930371065e13 / 0.510584971008e12 0.34774771016773e14 / 0.68928971086080e14; 0.18789143112277e14 / 0.45952647390720e14 -0.213895716727517e15 / 0.137857942172160e15 0.171024751153381e15 / 0.68928971086080e14 -0.8523669967037e13 / 0.2552924855040e13 0.485768751245399e15 / 0.137857942172160e15 -0.229158724354277e15 / 0.137857942172160e15; -0.51766014925489e14 / 0.689289710860800e15 0.202930494289627e15 / 0.689289710860800e15 -0.54332868549353e14 / 0.114881618476800e15 0.180479548146281e15 / 0.344644855430400e15 -0.1146942437956153e16 / 0.689289710860800e15 0.211001773091419e15 / 0.76587745651200e14;];



%     M=-(2*diag(ones(m-3,1),3)-27*diag(ones(m-2,1),2)+270*diag(ones(m-1,1),1)+270*diag(ones(m-1,1),-1)-27*diag(ones(m-2,1),-2)+2*diag(ones(m-3,1),-3)-490*diag(ones(m,1),0))/180;
    diags   = -3:3;
    stencil = -1/180*[2,-27,270,-490,270,-27,2];
    M = stripeMatrix(stencil, diags, m);

    M(1:6,1:6)=M_U;

    M(m-5:m,m-5:m)=rot90(  M_U ,2 );
    M=M/h;

    DS_U=[0.137e3 / 0.60e2 -5 5 -0.10e2 / 0.3e1 0.5e1 / 0.4e1 -0.1e1 / 0.5e1;];
    DS=sparse(m,m);
    DS(1,1:6)=DS_U;
    DS(m,m-5:m)=fliplr(DS_U);
    DS=DS/h;

    D2=H\(-M+DS);

%     d5=[-1 5 -10 10 -5 1];
%     t5=sum(abs(d5));
%     DD_5(1:1,1:6)=[d5];
%     DD_5(m:m,m-5:m)=[d5];
% 
%     % This works for wave eq.
%     % For studs interface in 1D no AD is needed.
%     ADD=7*h/(t5)*DD_5'*DD_5;

    e_1 = sparse(m,1);
    e_1(1)= 1;
    e_m = sparse(m,1);
    e_m(end)= 1;
    S_1 = -DS(1,:)';
    S_m =  DS(end,:)';

    Q = H*D1-(-(e_1*e_1') + (e_m*e_m'));
    M = -(H*D2-(-e_1*S_1' + e_m*S_m'));
end