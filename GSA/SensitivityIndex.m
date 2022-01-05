%To create the PRCC graph and boxplot:
clc
clear all
close all
%%% The sensitivity indices from the PRCC results, you'll need to change
%%% these if you run it on your own. The names on  line 8 correspond to the
%%% order of the variables in lines 11-15.
%'piM','miuM','betaM','piT','sigmaT','miuT','betaT','tauT','gammaT','miue','K', 'DM', 'DA', 'Pulses', 'tau'
% Columns 1 through 10

%     0.0481   -0.0060    0.0015    0.1583    0.5388   -0.0200    0.3750    0.6206    0.0145    0.0220
% 
%   Columns 11 through 15
% 
%     0.5721    0.0086    0.1054   -0.5324    0.9051
%    

%High Intensity
% sbetaM = 0.0015;
% sbetaT = 0.3750;
% smiuM = -0.0060;
% smiuT = -0.0200;
% sK = 0.5721;
% smiue = 0.0220;
% sgammaT = 0.0145;
% ssigmaT = 0.5388;
% stauT = 0.6206;
% spiT = 0.1583;
% spiM = 0.0481;
% sPulses = -0.5324;
% stau = 0.9051;
% sDM = 0.0086;
% sDA = 0.1054;


%'piM','miuM','betaM','piT','sigmaT','miuT','betaT','tauT','gammaT','miue','K', 'DM', 'DA', 'Pulses', 'tau'

%Columns 1 through 10
% 
%     0.0181    0.0532    0.0128    0.1662    0.6279   -0.0286    0.2675    0.6777    0.0362   -0.0197
% 
%   Columns 11 through 15
% 
%     0.6841   -0.0488    0.0452   -0.3660    0.9103
%Low Intensity
sbetaM = 0.0128;
sbetaT = 0.2675;
smiuM = 0.0532;
smiuT = -0.0286;
sK = 0.6841;
smiue = -0.0197;
sgammaT = 0.0362;
ssigmaT = 0.6279;
stauT = 0.6777;
spiT = 0.1662;
spiM = 0.0181;
sPulses = -0.3660;
stau = 0.9103;
sDM = -0.0488;
sDA = 0.0452;

SbetaM = 'Mice Transmission Probability';
SbetaT = 'Tick Transmission Probability';
SmiuM = 'Mice Death Rate';
SmiuT = 'Tick Death Rate';
SK = 'Carrying Capacity';
Smiue = 'Egg Inviability/Death Rate';
SgammaT = 'Adult Development Rate';
SsigmaT = 'Larvae Development Rate';
StauT = 'Nymph Development Rate';
SpiT = 'Tick Birth Rate';
SpiM = 'Mice Birth Rate';
SPulses = '# of Burns';
Stau = 'Time Between Burns';
SDM = 'Diffusion of Mice';
SDA = 'Diffusion of Adults';

sIN = [spiM smiuM sbetaM spiT ssigmaT smiuT sbetaT stauT sgammaT smiue sK sDM sDA sPulses stau];
%%% Definition of the parameters SmiuM Smue SgammaT StauT SpiT
IN = {SpiM SmiuM SbetaM SpiT SsigmaT SmiuT SbetaT StauT SgammaT Smiue SK SDM SDA SPulses Stau};
    

figure(2)
set(0,'DefaultAxesFontSize',20) 

barh(sIN);   %%% Plots the index on a bar graph
xlabel('Sensitivity index')
set(gca,'YLim',[1 40])       
    set(gca,'YTick',1:length(IN)) 
    set(gca,'YTickLabel',IN),  %%% Set the parameters definition on the Y-axis
    axis tight
colormap cool %summer HSV Summer 
title(['Low Intensity Burns'])
% ylabel('Description of Parameters')

%load matrix created from Ticks_LHS_PF for boxplots
LI_IN = cell2mat(struct2cell(load('In_LI.mat')));
HI_IN = cell2mat(struct2cell(load('In_HI.mat')));
LI_IN = LI_IN(:);
HI_IN = HI_IN(:);
isoutLI = isoutlier(LI_IN,'quartiles');
isoutHI = isoutlier(HI_IN,'quartiles');
LI_IN(isoutLI) = NaN;
HI_IN(isoutHI) = NaN;

figure(6)
set(gca, 'FontSize', 20)
boxplot([LI_IN,HI_IN],'Labels',{'Low Intensity','High Intensity'},'symbol', '')
ylim([0,2500])
% title('Dist. of Sum of Infectious Nymphs','fontsize',18)
