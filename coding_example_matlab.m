% This script construct a database based on the chosen coarse-grain model and assess it



clear all;
close all;

cm2eV = 1.2398e-4;
eV2j = 1.602e-19;
cm2j = cm2eV *eV2j;

kB = 1.38e-23;
NA = 6.022e23;

etrList = importdata('/hdd/tjpan/coarse_grain_model_in_DSMC/database/O2O_Varandas_RVT/Etr_list.input');
nEtr = size(etrList,1);
%A = importdata('/hdd/tjpan/coarse_grain_model_in_DSMC/database/O2O_Varandas_RVT/Erovib_oddOnly.input');
A = importdata('/hdd/tjpan/coarse_grain_model_in_DSMC/database/input/Erovib_sorted.input');
rvLevel = A(:,1:2);

nLevel = size(rvLevel,1);
ervList = [A(:, 1:2) A(:,4)];

ervList(:,3) = ervList(:,3)-A(1,6);

TList = 1000:1000:20000;
nTemp = size(TList,2);

nCoeff = 8;

folderInput =  '/hdd/tjpan/coarse_grain_model_in_DSMC/database/O2O_Varandas_RVT_organized/';
%folderInput =  '/hdd/tjpan/coarse_grain_model_in_DSMC/database/Test/';
folderOutput =  '/hdd/tjpan/coarse_grain_model_in_DSMC/database/O2O_Varandas_binned/';

%cd(folderInput);
%files =dir('crs*.*');


Idiss = find(ervList(:,1) ==47);% index of dissociation indicator

nQB = size(ervList,1);
IB = 1:1:Idiss-1;
IQ = Idiss:1:nQB;
%% convert i to g(group)

%%%%%%%%%%%%%%%%%%%%%
% set up the method of coarse-grain model
% uniform: fixed-size bins
% nth-order : variable-size bins
% tjpan : bi-variable-size bins

method = 'uniform';
%%%%%%%%%%%%%%%%%%%%%

% generate a rule to map state i to bin k
nB_k = 90;
nQ_k = 10;
nQB_k = nB_k +nQ_k;

Erange_B = ervList(Idiss,3)-0; %energy range of bound states
Erange_Q = ervList(IQ(end),3)-ervList(Idiss,3); %energy range of quasibound states

if strcmp(method,'uniform')
modelName = strcat(num2str(nB_k),':',num2str(nQ_k),'-eq');
elseif strcmp(method,'nth-order')
modelName = strcat(num2str(nB_k),':',num2str(nQ_k),'-var');
else
modelName = strcat(num2str(nB_k),':',num2str(nQ_k),'-tjpan');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
    case 'uniform'
        % constant energy size of bins
        deltaE_B = Erange_B/nB_k; 
        deltaE_Q = Erange_Q/nQ_k; 
        
        gridB = 0+ deltaE_B*(0:1:nB_k); 
        gridQ = ervList(Idiss,3)+deltaE_Q*(0:1:nQ_k);
    case 'nth-order'
        n =2;
        
        % constant energy size of bins
        deltaE_B = Erange_B/nB_k; 
        deltaE_Q = Erange_Q/nQ_k; 
        
        gridB = 0+ deltaE_B*(0:1:nB_k); 
        gridQ = ervList(Idiss,3)+deltaE_Q*(0:1:nQ_k);
        
        gridB = Erange_B*(gridB/Erange_B).^n;
        
        
     case 'tjpan'
        n =2;
        
        nB1_k = floor(nB_k/2);
        nB2_k = nB_k - nB1_k;
        Erange_B1 = Erange_B/2;
        Erange_B2 = Erange_B - Erange_B1;
        
        % constant energy size of bins
        deltaE_B = (Erange_B/2)/nB1_k; 
        deltaE_Q = Erange_Q/nQ_k; 
        
        idxB1 = 0:1:nB1_k;
        idxB2 = 1:1:nB2_k;
        gridQ = ervList(Idiss,3)+deltaE_Q*(0:1:nQ_k);
        
        gridB1 = Erange_B1*(idxB1/nB1_k).^n;
        gridB2 = Erange_B -(Erange_B2)*((nB2_k-idxB2)/nB2_k).^n;
        
        gridB = [gridB1 gridB2];
    otherwise
        error('required method is not defined')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i= 1 : nB_k
   index = find (ervList(:,3)>= gridB(i) & ervList(:,3)<= gridB(i+1)); 
   ervList(index,4) = i;
    
end

for i= 1 : nQ_k
   index = find (ervList(:,3)>= gridQ(i) & ervList(:,3)<= gridQ(i+1)); 
   ervList(index,4) = i + nB_k;
end

ervList(Idiss,4) = -1;

%% Bin energy


for i =1:nQB_k
    
    id = find (ervList(:,4)==i);
    gj = ervList(id,2)*2+1;%rotational degeneracy
    
    
ervBinList(i) = dot(ervList(id,3),gj)/sum(gj);
evibBinList(i) = dot(A(id,6),gj)/sum(gj);
erotBinList(i) = dot(A(id,7),gj)/sum(gj);
gjBin(i) = sum(gj);

end
%% visualization
x = 1:size(ervList,1);
figure;
for i =1 :nQB_k
    id = find (ervList(:,4)==i);
    if isempty(id)
        continue;
    else
    plot (x(id), ervList(id,3),'.');
    hold on
    line([x(max(id)) x(max(id))], [0 ervList(max(id),3)]);
    hold on
    line([0 x(max(id))], [ervList(max(id),3) ervList(max(id),3)]);
    hold on
    end
end

plot(zeros(nQB_k),ervBinList,'rx');
hold on
xlabel('i')
ylabel('E_{rv}[cm^{-1}]')

%% remove the bins containing no level
iBin =0;
idRemove =[];
for i =1 :nQB_k
    id = find (ervList(:,4)==i);
    if isempty(id)
        idRemove= [idRemove i];
    else
        iBin= iBin +1;
        ervList(id,4) = iBin;
        
    end
end
ervBinList(idRemove) =[];
evibBinList(idRemove)=[];
erotBinList(idRemove)=[];

nQB_k = iBin;

%% what is the equilibrium EDF of coarse-grain model

TList =[2000 5000 10000 20000];
legendList_diff = cell(length(TList),1);
outputbuf =[];

for j=1:length(TList)
T = TList(j);    
gjLevel = ervList(:,2)*2+1;%rotational degeneracy

%Z_torres = dot(gjBin,exp(-ErvBinList(:,3)*cm2j/k/T));
%Z_true = dot (gj,exp(-ErvList(:,3)*cm2j/k/T));

% plot EDF 
% _torres implies the results from coarse-grain equilibrium
% _tjpan implies the resultes from raw equilibrium, then grouped into bins

x_true = ervList(:,3)*cm2j;
y_true = gjLevel.*exp(-x_true/kB/T);
y_true_norm =y_true/sum(y_true);

x_torres = ervBinList*cm2j;
x_tjpan = x_torres;


for i= 1: nQB_k
    
  
   id = find (ervList(:,4)==i);
   gj = ervList(id,2)*2+1;%rotational degeneracy
   gk = sum(gj);
   
   y_torres(i) = gk*exp(-x_torres(i)/kB/T);
  
   x = ervList(id,3)*cm2j;
   
   y_tjpan(i) = sum(gj.*exp(-x/kB/T));
end
y_torres_norm =y_torres/sum(y_torres(~isnan(y_torres)));
y_tjpan_norm = y_tjpan/sum(y_tjpan);

% figure;
% semilogy (x_true,y_true/sum(y_true),'--',x_torres,y_torres_norm,'o',x_tjpan,y_tjpan_norm,'x');
% ylabel('n')
% 
% 
% figure;
% semilogy (x_torres,y_torres_norm,'o',x_tjpan,y_tjpan_norm,'x');
% ylabel('n')
% 

output = [y_torres_norm' y_tjpan_norm'];
rel_diff = max(y_torres_norm./y_tjpan_norm);


%


for i= 1: nQB_k
    
  
   id = find (ervList(:,4)==i);
   gj = ervList(id,2)*2+1;%rotational degeneracy
   gk(i) = sum(gj);
   
   y_torres(i) = exp(-x_torres(i)/kB/T);
  
   x = ervList(id,3)*cm2j;
   
   y_tjpan(i) = sum(gj.*exp(-x/kB/T))/gk(i);
end
y_torres_norm =y_torres/sum(y_torres);
y_tjpan_norm = y_tjpan/sum(y_tjpan);


figure;
semilogy (x_torres/eV2j,y_torres_norm,'o',x_tjpan/eV2j,y_tjpan_norm,'x');
legend('CG equilibrium','raw equilibrium')
xlabel('Erv[eV]')
ylabel('n/g')
title (strcat('EDF at T= ',num2str(T),' K'))
%legendList_EDF{i} = strcat('T= ',num2str(TList(i)),' K'); 

%relative difference
diff_torres = abs((y_torres_norm-y_tjpan_norm)./y_tjpan_norm*100);

% ratio
%diff_torres = abs((y_torres_norm)./y_tjpan_norm);


figure(3);
semilogy (x_torres/eV2j,diff_torres,'o-');

legendList_diff{j} = strcat('T= ',num2str(T),' K');
hold on


outputbuf = [outputbuf diff_torres']; 

diffArray(:,j) =diff_torres;



end
outputbuf = [x_torres'/eV2j outputbuf];

figure(3);
legend(legendList_diff,'Location','best');
xlabel('Erv[eV]')
ylabel('relative error [%]')
%ylabel('f_{CG}/f_{true}')
title (strcat('absolute relative error of EDF for',modelName) )
%ylim([0.01 10])

%% internal energy comparison

T = 100:10:20000;


% Gorden McBride database
R = 8.314510;

T_GM{1,1} = [200 1000];
T_GM{2,1} = [1000 6000];
T_GM{3,1} = [6000 20000];

a_GM{1,1} =[-3.425563420E+04 4.847000970E+02 1.119010961E+00 4.293889240E-03 -6.836300520E-07 -2.023372700E-09 1.039040018E-12];
b_GM{1,1} = [-3.391454870E+03 1.849699470E+01];

a_GM{2,1} = [-1.037939022E+06 2.344830282E+03 1.819732036E+00 1.267847582E-03 -2.188067988E-07...
       2.053719572E-11 -8.193467050E-16];
b_GM{2,1} = [-1.689010929E+04 1.738716506E+01];

a_GM{3,1} = [4.975294300E+08 -2.866106874E+05 6.690352250E+01 -6.169959020E-03 3.016396027E-07...
      -7.421416600E-12 7.278175770E-17];
   
b_GM{3,1} = [2.293554027E+06 -5.530621610E+02];

for i=1:length(T)
   if T(i)>=T_GM{1,1}(1,1) && T(i)<=T_GM{1,1}(1,2)
       flag=1;
   else if T(i)>T_GM{2,1}(1,1) && T(i)<=T_GM{2,1}(1,2)
        flag =2;
      
        else if T(i)>T_GM{3,1}(1,1) && T(i)<=T_GM{3,1}(1,2)
        flag =3;
            end
        end
   end
   a = a_GM{flag,1};
   b = b_GM{flag,1};
   
   t = T(i);
   Cp_GM(i)= R*(a(1)*t^-2+a(2)*t^-1+a(3)+a(4)*t+a(5)*t^2+a(6)*t^3+a(7)*t^4); % [J/mol-K]
   H_GM(i)= R*T(i)*(-a(1)*t^-2+a(2)*log(t)/t+a(3)+a(4)*t/2+a(5)*t^2/3+a(6)*t^3/4+a(7)*t^4/5+b(1)/t);%[J/mol] % w.r.t H(T= 298.15)=0
   S_GM (i)= R*(-a(1)*t^-2/2-a(2)*t^-1+a(3)*log(t)+a(4)*t+a(5)*t^2/2+a(6)*t^3/3+a(7)*t^4/4+b(2)); 
   U_GM(i) = H_GM(i) - R*t;
   Cv_GM(i)=Cp_GM(i) -R ;

end


%
Tref=298.15;
u_zero_true = sum(gjLevel.*x_true.*exp(-x_true/kB/Tref));

for i =1:length(T)
Z_true(i) = sum(gjLevel.*exp(-x_true/kB/T(i)));
u_true(i) = sum(gjLevel.*x_true.*exp(-x_true/kB/T(i)));%-u_zero_true;

Q = sum(gjLevel.*exp(-x_true/kB/T(i)));
A = sum(gjLevel.*x_true.*exp(-x_true/kB/T(i)));
B = sum(gjLevel.*x_true.^2.*exp(-x_true/kB/T(i)));

cv_true(i) = R*3/2+NA^2/R/T(i)^2*(B*Q-A^2)/Q^2; 
s_true(i) = R*(log(Q) + A/Q/kB/T(i));
end
u_true = u_true./Z_true*NA+3*R*T/2;


%
nk_eq_test = gk.*exp(-x_torres/kB/Tref);
u_zero_torres = sum(gk.*x_torres.*exp(-x_torres/kB/Tref));

for i =1:length(T)
Z_torres(i) = sum(gk.*exp(-x_torres/kB/T(i)));
u_torres(i) = sum(gk.*x_torres.*exp(-x_torres/kB/T(i)));%-u_zero_torres;

Q = sum(gk.*exp(-x_torres/kB/T(i)));
A = sum(gk.*x_torres.*exp(-x_torres/kB/T(i)));
B = sum(gk.*x_torres.^2.*exp(-x_torres/kB/T(i)));

cv_torres(i) = R*3/2+NA^2/R/T(i)^2*(B*Q-A^2)/Q^2; 
s_torres(i) = R*(log(Q) + A/Q/kB/T(i));
end
u_torres = u_torres./Z_torres*NA+3*R*T/2;

%
 for ii=1:nQB_k
      id = find (ervList(:,4)==ii);
      gj = ervList(id,2)*2+1;%rotational degeneracy
      nk_eq(ii) = sum(gj.*exp(-ervList(id,3)*cm2j/kB/Tref));    
 end
u_zero_tjpan = sum(nk_eq.*x_tjpan);

for i =1:length(T)

    for ii=1:nQB_k
      id = find (ervList(:,4)==ii);
      gj = ervList(id,2)*2+1;%rotational degeneracy
      nk_eq(ii) = sum(gj.*exp(-ervList(id,3)*cm2j/kB/T(i)));    
    end
Q = sum(nk_eq);
A = sum(nk_eq.*x_tjpan);
B = sum(nk_eq.*x_tjpan.^2);

Z_tjpan(i)=Q;
u_tjpan(i) = A;%-u_zero_tjpan;
cv_tjpan(i) = R*3/2+NA^2/R/T(i)^2*(B*Q-A^2)/Q^2; 

end
u_tjpan = u_tjpan./Z_tjpan*NA+3*R*T/2;
%
% figure;
% plot(T,H_GM/1000);
% xlabel('T[K]');
% ylabel('H[kJ/mol]');
% 
% figure;
% plot(T,Cp_GM);
% xlabel('T[K]');
% ylabel('Cp[J/mol-K]');
% 
% figure;
% plot(T,R*T*3/2/1000,'k--',T,R*T*5/2/1000,'k--',T,R*T*7/2/1000,'k--',T,U_GM/1000,'-',T,u_true/1000,'-');
% xlabel('T[K]');
% ylabel('U[kJ/mol]');
% legend('3/2RT','5/2RT','7/2RT','Gordon-McBride','state-specific','Location','best')

figure;
plot(T,R*T*5/2/1000,'k--',T,U_GM/1000,'--',T,u_true/1000,'--',T,u_torres/1000,'-.');
xlabel('T[K]');
ylabel('U[kJ/mol]');
legend('5/2RT','Gordon-McBride','state-specific','coarse-grain','Location','best')
title('specific internal energy')
%ylim([0 100])

figure;
plot(T,R*3/2*ones(length(T),1),'k--',T,R*5/2*ones(length(T),1),'k--',T,Cv_GM,'--',T,cv_true,'--',T,cv_torres,'-.');
xlabel('T[K]');
ylabel('Cv[J/mol-K]');
legend('3/2R','5/2R','Gordon-McBride','state-specific','coarse-grain','Location','best')
title('specific heat')
ylim([0 50])

figure;
plot(T,S_GM,'--',T,s_true,'--',T,s_torres,'-.');
xlabel('T[K]');
ylabel('S[J/mol-K]');
legend('Gordon-McBride','state-specific','coarse-grain','Location','best')
title('specific entropy')

%%
ans_GM = [T' U_GM' Cv_GM' S_GM'];
ans_true = [T' u_true' cv_true' s_true'];
ans_torres = [T' u_torres' cv_torres' s_torres'];
