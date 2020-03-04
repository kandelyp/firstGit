%% Prepare - Read state with Adiabatic projection 
%field parameters
clear all;
hbar=1./(2*pi);
Bz0=0.5;%%Tesla
Bz = [1  0 30 65].*0.001 + Bz0*6.15; % in GHz
sigmaBz=[1 1 1 1]*0.005;
U=[0 1]'; D=[1 0]';
Sigma{1}=[0 1; 1 0]; Sigma{2}=[0 -1j; 1j 0]; Sigma{3}=[1 0; 0 -1]; I=[1 0; 0 1];
Nspin=4;
Navg=64;
RtStab=[250 0 0].*1e-3;%evolve particular site adiabatically to prepare entangled state in GHz
RtStime=2000;%ns
dt=0.5;%ns
JaR=1; % way to measure entangled states with adread.
%initialize as S-T0 qubit
InitialState={'Qz', 'Qx'};%X: product state, Y: S+iT0 , Z = S
%measurement axis:
projAxes ={'X','Z','X'};%Z for diaread, X for adread and Y for projection along Y axis
%exchange Hamiltonian
Hex={};
for i=2:Nspin
    H=zeros(2^Nspin,2^Nspin);
    for j=1:3
        H = H+ kronProd({eye(2^(i-2),2^(i-2)), kronProd({Sigma{j},Sigma{j}}), eye(2^(Nspin-i),2^(Nspin-i))});
    end
    Hex{i-1}=0.25.*H;
end
H12 = 0.25*(kronProdOp('XXII') + kronProdOp('YYII') + kronProdOp('ZZII'));
for j=1:Navg
    
    Bz1=normrnd(Bz ,sigmaBz);
    %Hamiltonian for Zeeman interaction, assuming B=(0, 0, Bz);
    Hbz ={};HBz = zeros(2^Nspin,2^Nspin);
    for i=1:Nspin
        Hbz{i}= 0.5.*Bz1(i)*kronProd({eye(2^(i-1),2^(i-1)), Sigma{3}, eye(2^(Nspin-i),2^(Nspin-i))});
        HBz = HBz+Hbz{i};
    end
    
    
    %asign Spin at each site based on dBz between pair of spins
    Spin={};
    for i=2:2:Nspin
        if Bz1(i-1)>=Bz1(i)
            Spin{i-1}=U; Spin{i}=D;
            SpinGS(i-1)='U';SpinGS(i)='D';
        else
            Spin{i-1}=D; Spin{i}=U;
             SpinGS(i-1)='D';SpinGS(i)='U';
        end
    end
    
    %projection states:
    projMat={};
    for i=2:Nspin
        switch projAxes{i-1}
            case 'I'
                fprintf(['Identity operator is asigned as observable for ' num2str(i-1) 'th S-T0 qubit \n'])
            case 'X'
                temp = kronProd({Spin{i-1},Spin{i}});
            case 'Y'
                temp = kronProd({Spin{i-1},Spin{i}}) + 1j.*kronProd({Spin{i},Spin{i-1}});
            case 'Z'
                temp = kronProd({Spin{i-1},Spin{i}}) - kronProd({Spin{i},Spin{i-1}});
            otherwise
                error(['Wrong projection axis for ' num2str(i-1) ' S-T0 qubit\n']);
        end
       
        if projAxes{i-1}=='I'
            projMat{i-1}=eye(2^Nspin,2^Nspin);
        else
            %normalize temp
            temp = temp./norm(temp);
            projMat{i-1}= kronProd({eye(2^(i-2),2^(i-2)),temp*temp' , eye(2^(Nspin-i),2^(Nspin-i))});
        end
    end
    
    %Initialize wavefunction:
    psiInit=[];
    for i=2:2:Nspin
        psi=[];
        if strcmp(InitialState{i/2},'Qx') %product state
            psi = kronProd({Spin{i-1:i}});
        elseif strcmp(InitialState{i/2},'Qy') %S + iT0
            psi = kronProd({Spin{i-1:i}}) + 1j.*kronProd({Spin{i:-1:i-1}});
        elseif strcmp(InitialState{i/2},'Qz') %Singlet
            psi=kronProd({Spin{i-1:i}}) - kronProd({Spin{i:-1:i-1}});
        else
            error(['Wrong initial state for ' num2str(floor(i/2)) ' S-T0 qubit']);
        end
        if i==2
            psiInit=psi;
        else
            psiInit=kronProd({psiInit,psi});
        end
    end
    psiInit = psiInit./norm(psiInit);%normalize initial state
    
    %adiabatic ramp to prepare singlet in particular site:
    ind = find(RtStab);
    timeStep=0:dt:RtStime;%time steps for adiabatic ramp;
    N=length(timeStep);
    psi_t=psiInit;
    Jstep = linspace(0,RtStab(ind),N);
        if JaR ==1
            Jstep= flip(Jstep);
        end
    
    Prob={};
    for ii=1:N
        Uevol=expm(1j*(Hex{ind}.*Jstep(ii) + HBz).*dt./hbar);
        
        psi_t = Uevol*psi_t;
        
        %measure after every step:
        for i=2:Nspin
            Prob{i-1}(ii)=abs(psi_t'*projMat{i-1}*psi_t);
        end
    end
end
%% plot
figInd=10;figure(figInd);clf;
for i=1:length(Prob)
    plot(timeStep,Prob{i},'DisplayName',[ num2str(i) '    ' projAxes{i}]);hold on;
    lg=legend;
    title(lg,'Qubit  ProjAxis');
end
ylabel('P_S'); xlabel('evolution time(ns)');
sugercoatFig;


%% Adiabatic State Transfer with RampToSinglet 
%field parameters
figInd=1;
hbar=1./(2*pi);
Bz0=0.5;%%Tesla
Temp=0.1;%K
tm=5; T1=50; 
Bz = [ 6 0 10 16].*0.001+ 0*6.15; % in GHz
%Bz=[0 0 0 0; 10 0 50 0; 60 50 0 50; 10 0  50 100; 10 0 50 150;10  0 100 150; 110 100 50 0; 160 150 50 0; 160 150 100 0].*0.001+ 0.5*6.15;
sigmaBz=[1 1 1 1]*0.0025; 
U=[0 1]'; D=[1 0]';
Sigma{1}=[0 1; 1 0]; Sigma{2}=[0 -1j; 1j 0]; Sigma{3}=[1 0; 0 -1]; I=[1 0; 0 1];
S = (kron(U,D) - kron(D,U)).*0.5; T0=(kron(U,D) + kron(D,U)).*0.5; Tp= kron(U,U); Tm=kron(D,D);
loadFid =0.9; loadFid = [loadFid (1-loadFid)].^0.5; 
[tpL,~] = tripletFn(Bz0,Temp); tpL=tpL.^0.5;tpL=[1 0 0 0];
readFid=0.95;
swapFid=1;
%XX=0;
Q=[10 15 20];%Exchange Q factor
deltaJ=1./(sqrt(2)*pi*Q);%charge noise: dJ/J measured using exchange on 23

Nspin=4;
Navg=512;
nSteps =0:127; %ns %evolution time
dt=1;%integration time steps
piTime=5; %ns

jtab=[1 0  0; 0 1 0 ];
%jtab=eye(3,3);
rampFreq=120;
waitTime=0;%pulse collision wait time:
nonlinearJ=0;
alpha=1.0;

% Evolve to singlet or evolve from singlet using J ramp
RtStab=[0 0 0].*1e-3;%evolve particular site adiabatically to prepare entangled state in GHz
JaRJ=250e-3;%GHz Jump J for JaRdiaread
JaRtab=[];%JaRdiaread: [1-2 2-3 3-4], 1 for JaR, 0 otherwise, Jump and ramp for diaread
RtStime=[0.05 2000];%[ dt, Tevol] (ns), Ramp to singlet
%Ramp to singlet and Jump and Ramp use same time period and dt
%swap before read
swapRead=[3 2];%1 for 1-2, 2 for 2-3 , 3 for 3-4

%initialize as S-T0 qubit
InitialState={'Z', 'X'};%X: product state, Y: S+iT0 , Z = S, 'Tp for UU
%measurement axis:
projAxes ={'X','Z','Z'};%Z for diaread, X for adread and Y for projection along Y axis
%exchange Hamiltonian
Hex={};
for i=2:Nspin
    H=zeros(2^Nspin,2^Nspin);
    for j=1:3
        %         if XX &&j==3
        %             continue
        %         end
        
        H = H+ kronProd({eye(2^(i-2),2^(i-2)), kronProd({Sigma{j},Sigma{j}}), eye(2^(Nspin-i),2^(Nspin-i))});
    end
    Hex{i-1}=0.25.*H;
end
BStore=[];
tic
P_0=[]; P_1=[];P0=[]; P1=[];B=[];
for yp=1:length(rampFreq)
    %tic
    fprintf(['Running scan ' num2str(yp) ' of ' num2str(length(rampFreq)) '\n']);
    % jAQT=eye(3,3).*rampFreq(yp).*1e-3;%GHz
    jAQT=jtab.*rampFreq(yp).*1e-3;%GHz
    figure(55);clf;
    plot(jAQT.*1000);
    xticks([1:size(jAQT,1)]);
    %xticklabels({'1',num2str(nSteps)});
    legend({'J1','J2','J3'});
    title('AQT J ramp');xlabel('');
    ylabel('J (MHz)')
    pause(1);  
    
    
    % Then construct a ParforProgressbar object:
    ppm = ParforProgressbar(Navg);
    
    parfor j=1:Navg
        Bz1=normrnd(Bz ,sigmaBz);
        %Bz1=BStore(j,:);
        Prob0=[];Prob1=[];
        %Hamiltonian for Zeeman interaction, assuming B=(0, 0, Bz);
        Hbz ={};HBz = zeros(2^Nspin,2^Nspin);
        for i=1:Nspin
            Hbz{i}= 0.5.*Bz1(i).*kronProd({eye(2^(i-1),2^(i-1)), Sigma{3}, eye(2^(Nspin-i),2^(Nspin-i))});
            HBz = HBz+Hbz{i};
        end
        
        %asign Spin at each site based on dBz between pair of spins
        Spin={};SpinGS={};SpinErr={};
        for i=2:2:Nspin
            if Bz1(i-1)>Bz1(i)
                Spin{i-1}=U; Spin{i}=D;
                SpinErr{i-1} = D; SpinErr{i}=U;
                SpinGS{i-1}='U';SpinGS{i}='D';
            else
                Spin{i-1}=D; Spin{i}=U;
                SpinGS{i-1}='D';SpinGS{i}='U';
                SpinErr{i-1} = U; SpinErr{i}=D;
            end
            
        end
        
        %Initialize wavefunction:
        psiInit=[];
        for i=2:2:Nspin
            psi=[];
            if strcmp(InitialState{i/2},'X') %product state
                psi = loadFid(1).* kronProd({Spin{i-1:i}}) + loadFid(2).*(kronProd({SpinErr{i-1:i}}) +  Tm +Tp)./3;
                psi=psi./norm(psi);
            elseif strcmp(InitialState{i/2},'Y') %S + iT0
                psi = kronProd({Spin{i-1:i}}) + 1j.*kronProd({Spin{i:-1:i-1}});
                psi=psi./norm(psi);
            elseif strcmp(InitialState{i/2},'Z') %Singlet
                psi=loadFid(1).*S + loadFid(2).*(Tp + Tm + kronProd({SpinErr{i-1:i}}) + kronProd({Spin{i-1:i}}))./4;%kronProd({Spin{i-1:i}}) - kronProd({Spin{i:-1:i-1}});
                psi=psi./norm(psi);
            elseif strcmp(InitialState{i/2},'Tp')
                psi=tpL(1)*Tp +tpL(2)*S + tpL(3)*T0 + tpL(4)*Tm;
                psi=psi./norm(psi);
            else
                error(['Wrong initial state for ' num2str(floor(i/2)) ' S-T0 qubit']);
            end
            if i==2
                psiInit=psi;
                %psi=psi./norm(psi);
            else
                psiInit=kronProd({psiInit,psi});
            end
        end
        % psiInit = kronProd({S+1j*T0,D,U});
        
        psiInit = psiInit./norm(psiInit);%normalize initial state
        psi_t=psiInit;
        %[SpinGS{:}]
        %superpostion state at dot2
        %Spin{2}=(U + D)./sqrt(2); Spin{2}=Spin{2}./norm(Spin{2});
        %projection states:
        projMat={}; temp=[];
        for i=2:Nspin
            switch projAxes{i-1}
                case 'I'
                    fprintf(['Identity operator is asigned as observable for ' num2str(i-1) 'th S-T0 qubit \n'])
                case 'X'
                    temp = kronProd({Spin{i-1},Spin{i}});
                case 'Y' %works only for initialized paris i.e. 1-2 and 3-4, 2-3 might not always work this way
                    temp =0.5.* kronProd({Spin{i-1},Spin{i}}) + 1j.*kronProd({Spin{i},Spin{i-1}});
                case 'Z'
                    temp = S;%kronProd({Spin{i-1},Spin{i}}) - kronProd({Spin{i},Spin{i-1}});
                otherwise
                    error(['Wrong projection axis for ' num2str(i-1) ' S-T0 qubit\n']);
            end
            
            if projAxes{i-1}=='I'
                projMat{i-1}=eye(2^Nspin,2^Nspin);
            else
                %normalize temp
                temp = temp./norm(temp);
                projMat{i-1}= kronProd({eye(2^(i-2),2^(i-2)),temp*temp' , eye(2^(Nspin-i),2^(Nspin-i))});
            end
        end
        
        %            projMat{4}=kronProd({ Spin{3}*Spin{3}',eye(2^3,2^3)});
        %            projMat{5}=kronProd({eye(2,2), Spin{4}*Spin{4}',eye(2^2,2^2)});
        % %
        %projMat{4}=kronProd({ eye(2^3,2^3),U*U'});
        %projMat{5}=kronProd({eye(2^2,2^2), Spin{1}*Spin{1}',eye(2,2)});
        %         tmp=kronProd({Spin{2:4}});
        %         projMat{4}=kronProd({eye(2,2),tmp*tmp'});
        %         tmp = kronProd({Spin{4},Spin{3},Spin{2}});
        %         projMat{5}=kronProd({eye(2,2),tmp*tmp'});
        
        %adiabatic ramp to prepare singlet in particular site:
        ind = find(RtStab);
        if ind
            timeStep=0:RtStime(1):RtStime(2);%time steps for adiabatic ramp;
            N=length(timeStep);
            
            Jstep = linspace(0,RtStab(ind),N);
            Prob={};
            for ii=1:N
                Uevol=expm(1j*(Hex{ind}.*Jstep(ii) + HBz).*RtStime(1)./hbar);
                
                psi_t = Uevol*psi_t;
                
            end
        end
        
        psiInit=psi_t;
        if waitTime
            psiInit = expm(-1j.*HBz.*waitTime./hbar)*psiInit;
        end
        %adiabatic state transfer:
        for jj=1:length(nSteps)%
            wf1=psiInit; wf0=psiInit;
            if nSteps(jj)
                tEvol=dt:dt:nSteps(jj);
                AQTind = size(jAQT,1);
                for m=1:AQTind-1 %evolve from m th to m+1 th row of jAQT
                    J={};
                    for k=1:Nspin-1 %Js for each time evolution
                        if nonlinearJ
                            J{k}=rampFreq(yp).*1e-3.*linspace(jtab(m,k), jtab(m+1,k), length(tEvol)).^alpha;
                        else
                            J{k} = linspace(jAQT(m,k), jAQT(m+1,k), length(tEvol));
                            %                         J{k}= (jAQT(m,k)-jAQT(m+1,k)).*cos(pi.*(0:length(tEvol))./(2.*length(tEvol)));
                            %                         if jAQT(m,k)< jAQT(m+1,k)
                            %                             J{k}= abs(fliplr(J{k}));
                            %                         end
                        end
                        J{k}=normrnd(J{k}, J{k}.*deltaJ(k));
                    end
                    
                    
                    %time evolution:
                    for mm = 1:length(tEvol)
                        H1=zeros(2^Nspin,2^Nspin);
                        for k=1:Nspin-1
                            H1=H1+J{k}(mm).*Hex{k};
                        end
                        H0=H1;
                        H1=H1+HBz;
                        
                        U0=expm(-1j*H0.*dt./hbar);%evolution operator
                        U1=expm(-1j*H1.*dt./hbar);%evolution operator
                        wf0 = U0*wf0;
                        wf1 = U1*wf1;
                    end
                    
                end
            end
            
            %swap before read
            ind =swapRead;
            if ~ isempty(ind) && any(ind)
                for i=1:length(ind)
                    Hsr=zeros(size(HBz));
                    %for ii=1:length(ind)
                    Hsr = Hsr + Hex{ind(i)};
                    %end
                    Usr = expm(-1j.*waitTime.*HBz./hbar)* expm(1j.*(swapFid*pi*Hsr + (1.*HBz.*piTime./hbar)));
                    wf1 = Usr*wf1;
                    wf0 = Usr*wf0;
                end
            end
            
            %Jump and ramp for diaread
            ind = find(JaRtab);
            if ~isempty(ind)
                timeStep=0:RtStime(1):RtStime(2);%time steps for adiabatic ramp;
                N=length(timeStep);
                Jstep = linspace(0,JaRJ,N);
                Jstep= flip(Jstep);
                for ii=1:N
                    Hjar=zeros(size(HBz));
                    for iii=1:length(ind)
                        Hjar = Hjar + Jstep(ii).*Hex{ind(iii)};
                    end
                    Hjar = Hjar + HBz;
                    Ujar=expm(1j*Hjar.*RtStime(1)./hbar);
                    %psi_t = Uevol*psi_t;
                    wf1 = Ujar*wf1;
                    wf0=Ujar*wf0;
                end
            end
            
            if isempty(swapRead)
                wf1=expm(-1j.*HBz.*waitTime./hbar)*wf1;
            end
            
            for i=1:length(projMat) %for additional projection
                Prob0(i,jj)=abs(wf0'*projMat{i}*wf0);
                Prob1(i,jj)=abs(wf1'*projMat{i}*wf1);
            end
            
        end
        P0(j,:,:)=Prob0;
        P1(j,:,:)=Prob1;
        B(j,:)=Bz1;
        
        pause(100/Navg);
        % increment counter to track progress
        ppm.increment();
    end
    P_0(yp,:,:,:)=P0;
    P_1(yp,:,:,:)=P1;
    BStore(yp,:,:) = B;
    %
    %toc
end
% Delete the progress handle when the parfor loop is done (otherwise the timer that keeps updating the progress might not stop).
delete(ppm);
toc

 %% dBz
ln=1;
dBz=diff(BStore(:,:,:),[],3);

 %% 2D scan plot:
 % dBz=0 during evolution
Pavg=squeeze(nanmean(P_0,2));
figure(figInd);clf; pairs={'12','23','34','3','4'};projAxes{4}='S1'; projAxes{5}='S2';
for i=1:3
    ind=i;
    if i==2
        continue
    elseif i==3
        ind=2;
    end
    subplot(1,2,ind);
    imagesc(nSteps,rampFreq, squeeze(Pavg(:,i,:)));set(gca,'YDir','norm');colorbar;
    xlabel('RampTime (ns)'); ylabel('J (MHz)');title(['Pair-' pairs{i} ' ProjAxis-' projAxes{i}]);
    %caxis([0 1]);
end
sgtitle({[num2str(size(jAQT,1)+1) ' spin AQT'];...
    ['dBz=0 during J ramp'];...
    ['LoadType: ' [InitialState{:}] ' ' '     tm/T1= ' num2str(tm/T1)];...
    ['Load fidelity: S ' num2str(loadFid(1).^2) '  Tp ' num2str(round(tpL(1)^2,2))] })

%dBz during AQT
dBzdir=sign(dBz);
ind1=[];
for j=1:size(dBzdir,1)
ind1{j}=find(dBzdir(j,:,1)~=dBzdir(j,:,3));% .* dBzdir(:,2)~=dBzdir(:,1));
end

Pavg=[];
for j=1:length(ind1)
    Pavg(j,:,:)=squeeze(nanmean(squeeze(P_1(j,ind1{j},:,:)),1));
end
Pavg=squeeze(nanmean(P_1(:,:,:,:),2));
r=1-readFid; g=1-exp(-tm./T1);
%Pavg = (1-r).*Pavg +r.*(1-Pavg) +g.*(1-Pavg);
P1avg=Pavg;
figure(figInd+1);clf;
for i=1:3
    ind=i;
    if i==2
        continue
    elseif i==3
        ind=2;
    end
    subplot(1,2,ind);
    imagesc(nSteps,rampFreq, squeeze(Pavg(:,i,:)));set(gca,'YDir','norm');colorbar;
    xlabel('RampTime (ns)'); ylabel('J (MHz)');title(['Pair-' pairs{i} ' ProjAxis-' projAxes{i}]);
    %caxis([0 1]);
end
sgtitle({[num2str(size(jAQT,1)+1) ' spin AQT'];...
    [sprintf('Bz: 3.07GHz+ [%3.1f %3.1f %3.1f %3.1f ] MHz ', (Bz-3.0750).*1e3) ' ' sprintf('SigmaBz:  [%3.1f %3.1f %3.1f %3.1f ] MHz ', 1e3.*sigmaBz(:)) ];...
    ['LoadType: ' [InitialState{:}] ' ' '     tm/T1= ' num2str(tm/T1) 'Load fidelity: S ' num2str(loadFid(1).^2) '  Tp ' num2str(round(tpL(1)^2,2)) 'Readout fidelity: ' num2str(readFid) ' swapFid : ' num2str(swapFid) ' waitTime: ' num2str(waitTime) 'ns'] })
sugercoatFig;

Pavg=squeeze(nanmean(P_1,2));

figure(figInd+4);clf;
for i=4:5
    subplot(1,2,i-3);
    imagesc(nSteps,rampFreq, (1-(squeeze(Pavg(:,i,:)))));set(gca,'YDir','norm');colorbar;
    xlabel('RampTime (ns)'); ylabel('J (MHz)');title(['Pair-' pairs{i} ' ProjAxis-' projAxes{i}]);
    %caxis([0 1]);
    set(gca,'ColorScale','log');
end
sgtitle({[num2str(size(jAQT,1)+1) ' spin AQT'];...
    [sprintf('Bz: 3.07GHz+ [%3.1f %3.1f %3.1f %3.1f ] MHz ', (Bz-3.0750).*1e3) ' ' sprintf('SigmaBz:  [%3.1f %3.1f %3.1f %3.1f ] MHz ', 1e3.*sigmaBz(:)) ];...
    ['LoadType: ' [InitialState{:}] ' ' '     tm/T1= ' num2str(tm/T1)];...
    ['Load fidelity: S ' num2str(loadFid(1).^2) '  Tp ' num2str(round(tpL(1)^2,2))] })
%% linecut
ind=find(rampFreq==100);
p=squeeze(P1avg(ind,[1 3],:));
figure(figInd+3);clf;
plot(squeeze(p'));hold on;

% Data and Simulation comparison
figure(556);clf;
xval = (0:size(p,2)-1);
plot(xval,p','--'); hold on;
plot(xval, lnctData(:,1:size(p,2))');
legend({'A','B','A','B'});
xlabel('ramp time (ns)'); ylabel('P_S'); 
grid on; sugercoatFig;

Bzs=(Bz-3.0750).*1e3;
title({'Solid line: Experiment, Dashed line: Simulation , J=100 MHz';...
    [sprintf('Bz: 3.07GHz+ [%3.1f %3.1f %3.1f %3.1f ] MHz ', Bzs(:)) ' ' sprintf('SigmaBz:  [%3.1f %3.1f %3.1f %3.1f ] MHz ', 1e3.*sigmaBz(:)) ];...
        ['LoadType: ' [InitialState{:}] ' ' '     tm/T1= ' num2str(tm/T1)    'Load fidelity: S ' num2str(loadFid(1).^2) '  Tp ' num2str(round(tpL(1)^2,2)) ' SwapFid ' num2str(swapFid)] })


%% 1D scan plot
dBz=diff(squeeze(BStore),[],2);
dBzdir=sign(dBz);
ind=find(dBzdir(:,1)~=dBzdir(:,3));% .* dBzdir(:,2)~=dBzdir(:,1));
projAxes{4}='234';projAxes{5}='432';
p=[];
figInd=6;
figure(figInd);hold on;
leg={'1-2','2-3','3-4','3','4'};
j=1;
for i=1:size(P1,2)
    if any([ 2 4 5 ]==i)
        continue
    end
    %p(j,:)= squeeze(nanmean(P1(ind(:),i,:),1));
    p(j,:)= squeeze(nanmean(P1(ind,i,:),1));
    r=1-readFid; g=1-exp(-tm./T1);
    p(j,:) = (1-r).*p(j,:) +r.*(1-p(j,:)) +g.*(1-p(j,:));
    plot(nSteps,p(j,:) ,'DisplayName',[leg{i} '         ' projAxes{i}]);
    ylabel('P_S');xlabel('Ramp time (ns)');
    j=j+1;
end

line([0 nSteps(end)],[0.5 0.5],'DisplayName','0.5','color','g');
title('with dBz');
ylim([0 1]);
legend show;
lg = legend;
title(lg,[' Spin    ProjAxis']);
%
figure(figInd+1);clf;hold on;
p=squeeze(nanmean(P0(:,:,:),1));
for i=1:size(P0,2)
    if any([2 4 5]==i)
        continue
    end
   % plot(nSteps, squeeze(nanmean(P0(ind,i,:),1)),'DisplayName',[leg{i} '         ' projAxes{i}] );
     plot(nSteps, squeeze(nanmean(P0(:,i,:),1)),'DisplayName',[leg{i} '         ' projAxes{i}] );
end
line([0 nSteps(end)],[0.5 0.5],'DisplayName','0.5','color','g');
title('dBz=0 during Jramp');
legend show;
lg = legend;
title(lg,['       Spin    ProjAxis']);
%sugercoatFig;
Bzs=(Bz-3.0750).*1e3;
sgtitle({'AQT';...
    [sprintf('Bz: 3.07GHz+ [%3.1f %3.1f %3.1f %3.1f ] MHz ', Bzs(:)) ' ' sprintf('SigmaBz:  [%3.1f %3.1f %3.1f %3.1f ] MHz ', 1e3.*sigmaBz(:)) ];...
        ['LoadType: ' [InitialState{:}] ' ' '     tm/T1= ' num2str(tm/T1)];...
        ['Load fidelity: S ' num2str(loadFid(1).^2) '  Tp ' num2str(round(tpL(1)^2,2)) ' SwapFid' num2str(swapFid)] })

% figure(figInd+1);clf;
% plot(jAQT.*1000);
% xticks([1:size(jAQT,1)]);
% %xticklabels({'1',num2str(nSteps)});
% legend({'J1','J2','J3'});
% title('AQT J ramp');xlabel('');
% ylabel('J (MHz)');
%sugercoatFig;
%% Data and Simulation comparison
figure(556);clf;
xval = (0:size(p,2)-1);
plot(xval,p','--'); hold on;
plot(xval, lnctData(:,1:size(p,2))');
grid on; %sugercoatFig;

Bzs=(Bz-3.0750).*1e3;
title({'Solid line: Experiment, Dashed line: Simulation , J=100 MHz';...
    [sprintf('Bz: 3.07GHz+ [%3.1f %3.1f %3.1f %3.1f ] MHz ', Bzs(:)) ' ' sprintf('SigmaBz:  [%3.1f %3.1f %3.1f %3.1f ] MHz ', 1e3.*sigmaBz(:)) ];...
        ['LoadType: ' [InitialState{:}] ' ' '     tm/T1= ' num2str(tm/T1)    'Load fidelity: S ' num2str(loadFid(1).^2) '  Tp ' num2str(round(tpL(1)^2,2)) ' SwapFid ' num2str(swapFid)] })
   
%% Plot dBz effect test
h1=figure(111);clf; 
ps=squeeze(P_1(:,:,[1 3],:));
clr = colormap(parula(9));
clr={'r','g','b','r','g','b','r','g','b'};
Bz=[0 0 0 0; 10 0 50 0; 60 50 0 50; 10 0  50 100; 10 0 50 150;10  0 100 150; 110 100 50 0; 160 150 50 0; 160 150 100 0].*0.001+ 0.5*6.15;
for i=1:9
    p=squeeze(ps(i,:,:));
     r=1-readFid; g=1-exp(-tm./T1);
    p = (1-r).*p +r.*(1-p) +g.*(1-p);
    k=ceil(i/3);   
    subplot(1,3,k);hold on;
    a(i)=plot(nSteps,p(1,:),'DisplayName',sprintf('%0.0f    %0.0f   %0.0f ',diff(Bz(i,1:4))*1e3),'Color',clr{i}); hold on;
    plot(nSteps, p(2,:),':','DisplayName',sprintf('%0.0f    %0.0f   %0.0f ',diff(Bz(i,1:4))*1e3),'Color',clr{i}); hold on;
    %legend show;
    xlabel('t (ns)'); ylabel('P_S');
    
end
for i=1:3
    subplot(1,3,i);
   legend(a((i-1)*3+1:i*3));
end
sugercoatFig;


sgtitle({[' Solid line:A(Z)  Dotted Line: B(X)   LoadType: ' [InitialState{:}] ' ' '     tm/T1= ' num2str(tm/T1) '   Load fidelity: S ' num2str(loadFid(1).^2) '  Tp ' num2str(round(tpL(1)^2,2)) ' MeasFid: ' num2str(readFid) ] })
width=12; height=0.33*width;
set(h1,'Units','inches','Position',[10 5 width height])
figure(112);clf;
plot(0:1,jAQT.*1e3);
xticks([0 1]); xticklabels({'0','127'});
xlabel('t(ns)'); ylabel('maxJ');
legend({'J1','J2','J3'}); sugercoatFig;
    
%% Resonant peaks and Infidelity 
p= squeeze(P_0);
%p=squeeze(P0(1,4,:));
J=rampFreq*1e-3;
xval=nSteps*J*hbar;
figure(100);clf; hold on;
r=1-readFid; g=1-exp(-tm./T1);
    
p= squeeze(P_0(1,1,1,:));
%p= (1-r).*p +r.*(1-p) +g.*(1-p);
plot(xval,abs(p),'--','DisplayName','dBz=0');hold on;

p= squeeze(P_1(1,1,1,:));
%p= (1-r).*p +r.*(1-p) +g.*(1-p);
plot(xval,abs(p),':','DisplayName','Bz= [20 0 60 90]+B0]');hold on;
legend show; xlabel('JT $\hbar$','Interpreter', 'latex'); ylabel('F');
 set(gca,'YScale','log');
% [pk, loc]=findpeaks(abs(log(1-p(4,:))),xval);
% 
% figure(4);clf;plot(diff(loc))
%% save data


%% Evolve- AST - Diaread on both sides:
%field parameters
clear all;
hbar=1./(2*pi);
Bz0=0.5;%%Tesla
Temp = 0.05;%temperature, K
Bz = [78 0 60 146].*0.001 + Bz0*6.15; % in GHz
sigmaBz=[1 1 1 1]*0.012;
U=[0 1]'; D=[1 0]';
Sigma{1}=[0 1; 1 0]; Sigma{2}=[0 -1j; 1j 0]; Sigma{3}=[1 0; 0 -1]; I=[1 0; 0 1];
S = (kron(U,D) - kron(D,U)).*0.5; T0=(kron(U,D) + kron(D,U)).*0.5; Tp= kron(U,U); Tm=kron(D,D);

%load error
loadFid = 0.9; loadFid = [loadFid 1-loadFid].^0.5;
[tpL,~] = tripletFn(Bz0,Temp); tpL=tpL.^0.5;

Nspin=4;
Navg=1;
nSteps = 127; %ns %evolution time
dt=1;%integration time steps
piTime=5; %ns
%AST ramp J: from ith row to i+1 th row
jAQT=eye(3,3).*0.1;%GHz


% Evolve to singlet or evolve from singlet using J ramp
RtStab=[0 0 0].*1e-3;%evolve particular site adiabatically to prepare entangled state in GHz
JaRJ=250e-3;%GHz Jump J for JaRdiaread
JaRtab=[0 0 0];%JaRdiaread: [1-2 2-3 3-4], 1 for JaR, 0 otherwise, Jump and ramp for diaread
RtStime=[0.05 1000];%[ dt, Tevol] (ns), Ramp to singlet
%Ramp to singlet and Jump and Ramp use same time period and dt

%swap before read
swapRead=[];%1 for 1-2, 2 for 2-3 , 3 for 3-4

%evolve state before AST
evolTab = [1, 127, 2];% [control (1 or 0), nsteps, dt(ns)]
%initialize as S-T0 qubit
InitialState={'Z', 'X'};%X: product state, Y: S+iT0 , Z = S, 'Tp for UU
%measurement axis:
projAxes ={'X','Z','Z'};%Z for diaread, X for adread and Y for projection along Y axis
%exchange Hamiltonian
Hex={};

for i=2:Nspin
    H=zeros(2^Nspin,2^Nspin);
    for j=1:3
        H = H+ kronProd({eye(2^(i-2),2^(i-2)), kronProd({Sigma{j},Sigma{j}}), eye(2^(Nspin-i),2^(Nspin-i))});
    end
    Hex{i-1}=0.25.*H;
end

for j=1:Navg
    Prob0=[];Prob1=[];
    temp=[];
    Bz1=normrnd(Bz ,sigmaBz);
    %Hamiltonian for Zeeman interaction, assuming B=(0, 0, Bz);
    Hbz ={};HBz = zeros(2^Nspin,2^Nspin);
    for i=1:Nspin
        Hbz{i}= 0.5.*Bz1(i)*kronProd({eye(2^(i-1),2^(i-1)), Sigma{3}, eye(2^(Nspin-i),2^(Nspin-i))});
        HBz = HBz+Hbz{i};
    end
    
    %asign Spin at each site based on dBz between pair of spins
    Spin={};SpinGS={};SpinErr={};
    for i=2:2:Nspin
        
        if Bz1(i-1)>=Bz1(i)
            Spin{i-1}=U; Spin{i}=D;
            SpinErr{i-1}=D;  SpinErr{i}=U;
            SpinGS{i-1}='U';SpinGS{i}='D';
        else
            Spin{i-1}=D; Spin{i}=U;
            SpinErr{i-1}=U;  SpinErr{i}=D;
            SpinGS{i-1}='D';SpinGS{i}='U';
        end
        
    end
    %projection states:
    projMat={};
    for i=2:Nspin
        switch projAxes{i-1}
            case 'I'
                fprintf(['Identity operator is asigned as observable for ' num2str(i-1) 'th S-T0 qubit \n'])
            case 'X'
                temp = kronProd({Spin{i-1},Spin{i}});
            case 'Y'
                temp = kronProd({Spin{i-1},Spin{i}}) + 1j.*kronProd({Spin{i},Spin{i-1}});
            case 'Z'
                temp = S;%kronProd({Spin{i-1},Spin{i}}) - kronProd({Spin{i},Spin{i-1}});
            otherwise
                error(['Wrong projection axis for ' num2str(i-1) ' S-T0 qubit\n']);
        end
        
        if projAxes{i-1}=='I'
            projMat{i-1}=eye(2^Nspin,2^Nspin);
        else
            %normalize temp
            temp = temp./norm(temp);
            projMat{i-1}= kronProd({eye(2^(i-2),2^(i-2)),temp*temp' , eye(2^(Nspin-i),2^(Nspin-i))});
        end
    end
    
    %Initialize wavefunction:
    psiInit=[];
    for i=2:2:Nspin
        psi=[];
        if strcmp(InitialState{i/2},'X') %product state
            psi = loadFid(1).* kronProd({Spin{i-1:i}}) + loadFid(2).*(kronProd({SpinErr{i-1:i}}) + Tm +Tp)/3;
            psi=psi./norm(psi);
        elseif strcmp(InitialState{i/2},'Y') %S + iT0
            psi = kronProd({Spin{i-1:i}}) + 1j.*kronProd({Spin{i:-1:i-1}});
            psi=psi./norm(psi);
        elseif strcmp(InitialState{i/2},'Z') %Singlet
            psi=loadFid(1).*S + loadFid(2).*(Tp + Tm + T0)./3;%kronProd({Spin{i-1:i}}) - kronProd({Spin{i:-1:i-1}});
            psi=psi./norm(psi);
        elseif strcmp(InitialState{i/2},'Tp')
            psi=tpL(1)*Tp +tpL(2)*S + tpL(3)*T0 + tpL(4)*Tm;
            psi=psi./norm(psi);
        else
            error(['Wrong initial state for ' num2str(floor(i/2)) ' S-T0 qubit']);
        end
        if i==2
            psiInit=psi;
            psi=psi./norm(psi);
        else
            psiInit=kronProd({psiInit,psi});
        end
    end
    psiInit = psiInit./norm(psiInit);%normalize initial state
    psi_t=psiInit;
    
    %adiabatic ramp to prepare singlet in particular site:
    ind = find(RtStab);
    if ind
        timeStep=0:RtStime(1):RtStime(2);%time steps for adiabatic ramp;
        N=length(timeStep);
        
        Jstep = linspace(0,RtStab(ind),N);
        Prob={};
        for ii=1:N
            Uevol=expm(1j*(Hex{ind}.*Jstep(ii) + HBz).*RtStime(1)./hbar);
            
            psi_t = Uevol*psi_t;
            
            %measure after every step:
            %         for i=2:Nspin
            %             Prob{i-1}(ii)=abs(psi_t'*projMat{i-1}*psi_t);
            %         end
        end
    end
    
    tEvol=0:dt:nSteps;AQTind = size(jAQT,1);
    for jj=1:evolTab(2)%
        psiInit=psi_t;
        %evolve
        if evolTab(1)
            psiInit =  expm(1j*HBz*(jj-1)*evolTab(3)./hbar)*psiInit;
        end        
        %AST        
        
        wf1=psiInit; wf0=psiInit;
        for m=1:AQTind-1 %evolve from m th to m+1 th row of jAQT
            J={};
            for k=1:Nspin-1 %Js for each time evolution
                J{k} = linspace(jAQT(m,k), jAQT(m+1,k), length(tEvol));
            end
            %time evolution:
            for mm = 1:length(tEvol)                
                H1=zeros(2^Nspin,2^Nspin);
                for k=1:Nspin-1
                    H1=H1+J{k}(mm).*Hex{k};
                end
                H0=H1;
                H1=H1+HBz;
                
                U0=expm(1j*H0.*dt./hbar);%evolution operator
                U1=expm(1j*H1.*dt./hbar);%evolution operator
                wf0 = U0*wf0;
                wf1 = U1*wf1;
            end
        end       

        %swap before read
        ind =swapRead;
        if ~ isempty(ind)
            Hsr=zeros(size(HBz));
            for ii=1:length(ind)
                Hsr = Hsr + Hex{ind(ii)};
            end
            Usr =  expm(1j.*(pi*Hsr + (HBz.*piTime./hbar)));
            wf1 = Usr*wf1;
            wf0 = Usr*wf0;
        end
        
        %Jump and ramp for diaread
        ind = find(JaRtab);
        if ~isempty(ind)
            timeStep=0:RtStime(1):RtStime(2);%time steps for adiabatic ramp;
            N=length(timeStep);
            Jstep = linspace(0,JaRJ,N);
            Jstep= flip(Jstep);
            for ii=1:N
                Hjar=zeros(size(HBz));
                for iii=1:length(ind)
                    Hjar = Hjar + Jstep(ii).*Hex{ind(iii)};
                end
                Hjar = Hjar + HBz;
                Ujar=expm(1j*Hjar.*RtStime(1)./hbar);
                %psi_t = Uevol*psi_t;
                wf1 = Ujar*wf1;
                wf0=Ujar*wf0;
            end
        end
        
        %measure state of all qubits
        for i=1:Nspin-1
            Prob0(i,jj)=abs(wf0'*projMat{i}*wf0);
            Prob1(i,jj)=abs(wf1'*projMat{i}*wf1);
        end       
    end
    Bstore(j,:)=Bz1;
    P0(j,:,:)=Prob0;
    P1(j,:,:)=Prob1;
    
end

%% plot

tm=4;T1=50;
readFid=0.95;

figInd=1;
figure(figInd);clf;
leg={'1-2','2-3','3-4'};
subplot(1,2,1);hold on;

j=1;
for i=1:size(P1,2)
    if any([ 2 4 5 ]==i)
        continue
    end
    %p(j,:)= squeeze(nanmean(P1(ind(:),i,:),1));
    p(j,:)= squeeze(nanmean(P1(:,i,:),1));
    r=1-readFid; g=1-exp(-tm./T1);
    p(j,:) = (1-r).*p(j,:) +r.*(1-p(j,:)) +g.*(1-p(j,:));
    plot((0:evolTab(2)-1).*evolTab(3),p(j,:) ,'DisplayName',[leg{i} '         ' projAxes{i}]);
   % ylabel('P_S');xlabel('Ramp time (ns)');
    j=j+1;
end

line([0 nSteps(end)],[0.5 0.5],'DisplayName','0.5','color','g');

line([0 evolTab(2)*evolTab(3)],[0.5 0.5],'DisplayName','0.5','color','g');
title('dBz!=0 during AQT J ramp');ylabel('P_S'); xlabel('dBz evolution time(ns)');
lg=legend;
ylim([0 1]);title(lg,'       Spin    ProjAxis');
%xlim([0 30]);
subplot(1,2,2);hold on;
for i=1:size(P0,2)
    
plot((0:evolTab(2)-1).*evolTab(3), squeeze(nanmean(P0(:,i,:),1)),'DisplayName',[leg{i} '         ' projAxes{i}] );
end
line([0 evolTab(2)*evolTab(3)],[0.5 0.5],'DisplayName','0.5','color','g');
title('dBz=0 during AQT J ramp');
legend show;
lg = legend;
title(lg,'       Spin    ProjAxis');
sugercoatFig;
Bzs=(Bz-3.0750).*1e3;
sgtitle({'Evolve around dBz and AQT';...
    [sprintf('Bz: 3.07GHz+ [%3.1f %3.1f %3.1f %3.1f ] MHz ', Bzs(:)) ' ' sprintf('SigmaBz:  [%3.1f %3.1f %3.1f %3.1f ] MHz ', 1e3.*sigmaBz(:)) ' MHz'];...
        ['LoadType: ' [InitialState{:}] ' ' '     tm/T1= ' num2str(tm/T1)];...
        ['Load fidelity: S ' num2str(loadFid(1).^2) '  Tp ' num2str(round(tpL(1)^2,2))] })
figure(figInd+1);clf;

plot(jAQT.*1000);
xticks([1:size(jAQT,1)]);
xticklabels({'1',num2str(nSteps)});
legend({'J1','J2','J3'});
title('AQT J ramp');xlabel(' Time (ns)');
ylabel('J (MHz)');
sugercoatFig;

%% 2D plots
figure(figInd+2);clf; 
frq=[];
for i=1:64
    B(i,:)=Bstore{i};
end
dBzs=abs(diff(B,[], 2));

sides={'A','B'};
for i = 1:2
    subplot(2,2,i);
    ind=i;
    if i>1
        ind=3;
    end
    p=squeeze(P1(:,ind,:));
    p=1- ((1-p).*exp(-tm/T1));
    imagesc(p);colorbar;title(['Data- Side' sides{i}]);
    set(gca,'YDir','norm'); xlabel('dBz evolution time(ns)');ylabel('repetition');
    p=p-nanmean(p,1);
    f=abs(fft(p,[],2));  
    L=size(p,2);
    xval=1e3./evolTab(3) .* (0:L/2-1)./L;
    subplot(2,2,2+i);
    imagesc(xval,1:size(f,1),  f(:,1:length(xval)));colorbar;
    %imagesc(f);colorbar;
    title(['FFT- Side' sides{i}]);
    set(gca,'YDir','norm'); xlabel('Frequency (MHz)');    ylabel('repetition');    
    [~,ind]= max(f(:,1:L/2),[],2);
    frq(i,:)=xval(ind); 
    hold on;
    plot(frq(i,:)',1:size(f,1), 'or');
end

figure(figInd+3);clf; hold on;
plot(1e3.*dBzs(:,[1 3]),'o'); 
plot(frq'); 
xlabel('repetition'); ylabel('Frequency (MHz)');
legend({'dBz12','dBz34','Side A','SideB'});sugercoatFig;