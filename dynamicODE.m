%% modified fibroblast model ODE
function [dydt]= dynamicODE(t,y,params) 
% Assign names for parameters 
[rpar,tau,ymax,speciesNames,KI,Inputs,inptIns,hrs]=params{:}; 

w = rpar(1,:);
n = rpar(2,:);
EC50 = rpar(3,:);
dose = rpar(4,:);
drugBinding = rpar(5,:); 
drugAgonism = rpar(6,:); 
w(inptIns) = interp1(hrs,Inputs',t);

w(inptIns)


rpar = [w;n;EC50;dose;drugBinding;drugAgonism];

AngII = 1; 
inputAngII = 2; 
AngIImRNA = 3; 
AT1R = 4; 
AGT = 5; 
ACE = 6; 
NOX = 7; 
ROS = 8; 
ET1 = 9; 
inputET1 = 10; 
ET1mRNA = 11; 
ETAR = 12; 
DAG = 13; 
PKC = 14; 
TRPC = 15; 
NE = 16; 
inputNE = 17; 
BAR = 18; 
Forskolin = 19; 
inputForskolin = 20; 
AC = 21; 
cAMP = 22; 
PKA = 23; 
CREB = 24; 
CBP = 25; 
TGFB = 26; 
inputTGFB = 27; 
TGFB1R = 28; 
smad3 = 29; 
smad7 = 30; 
TGFBmRNA = 31; 
BAMBI = 32; 
PDGF = 33; 
inputPDGF = 34; 
PDGFR = 35; 
NP = 36; 
inputNP = 37; 
NPRA = 38; 
cGMP = 39; 
PKG = 40; 
mechanical = 41; 
inputMechanical = 42; 
B1int = 43; 
Rho = 44; 
ROCK = 45; 
Ca = 46; 
calcineurin = 47; 
NFAT = 48; 
IL6 = 49; 
inputIL6 = 50; 
IL6mRNA = 51; 
gp130 = 52; 
STAT = 53; 
IL1 = 54; 
inputIL1 = 55; 
IL1RI = 56; 
TNFa = 57; 
inputTNFa = 58; 
TNFaR = 59; 
NFKB = 60; 
PI3K = 61; 
Akt = 62; 
p38 = 63; 
TRAF = 64; 
ASK1 = 65; 
MKK3 = 66; 
PP1 = 67; 
JNK = 68; 
abl = 69; 
Rac1 = 70; 
MEKK1 = 71; 
MKK4 = 72; 
ERK = 73; 
Ras = 74; 
Raf = 75; 
MEK1 = 76; 
FAK = 77; 
epac = 78; 
Factin = 79; 
FA = 80; 
migration = 81; 
cmyc = 82; 
CTGF = 83; 
proliferation = 84; 
SRF = 85; 
EDAFN = 86; 
aSMA = 87; 
AP1 = 88; 
TIMP1 = 89; 
TIMP2 = 90; 
PAI1 = 91; 
proMMP14 = 92; 
proMMP1 = 93; 
proMMP2 = 94; 
proMMP9 = 95; 
MMP1 = 96; 
MMP2 = 97; 
MMP9 = 98; 
MMP14 = 99; 
periostin = 100; 
CImRNA = 101; 
CIIImRNA = 102; 
CI = 103; 
CIII = 104; 
LOX = 105; 
xlinkFibers = 106; 
contraction = 107; 
dydt = zeros(107,1); 
dydt(AngII) = (OR(act(y(inputAngII),rpar(:,12)),act(y(AngIImRNA),rpar(:,27)))*ymax(AngII) - y(AngII))/tau(AngII); 
dydt(inputAngII) = (rpar(1,1)*ymax(inputAngII) - y(inputAngII))/tau(inputAngII); 
dydt(AngIImRNA) = (AND(rpar(:,26),act(y(AGT),rpar(:,26)),act(y(ACE),rpar(:,26)))*ymax(AngIImRNA) - y(AngIImRNA))/tau(AngIImRNA); 
dydt(AT1R) = (act(y(AngII),rpar(:,34))*ymax(AT1R) - y(AT1R))/tau(AT1R); 
dydt(AGT) = (AND(rpar(:,46),inhib(y(AT1R),rpar(:,46)),act(y(p38),rpar(:,46)),inhib(y(JNK),rpar(:,46)))*ymax(AGT) - y(AGT))/tau(AGT); 
dydt(ACE) = (act(y(TGFB1R),rpar(:,69))*ymax(ACE) - y(ACE))/tau(ACE); 
dydt(NOX) = (OR(act(y(AT1R),rpar(:,35)),act(y(TGFB1R),rpar(:,132)))*ymax(NOX) - y(NOX))/tau(NOX); 
dydt(ROS) = (OR(act(y(NOX),rpar(:,36)),act(y(ETAR),rpar(:,56)))*ymax(ROS) - y(ROS))/tau(ROS); 
dydt(ET1) = (OR(act(y(inputET1),rpar(:,20)),act(y(ET1mRNA),rpar(:,33)))*ymax(ET1) - y(ET1))/tau(ET1); 
dydt(inputET1) = (rpar(1,9)*ymax(inputET1) - y(inputET1))/tau(inputET1); 
dydt(ET1mRNA) = (act(y(AP1),rpar(:,32))*ymax(ET1mRNA) - y(ET1mRNA))/tau(ET1mRNA); 
dydt(ETAR) = (act(y(ET1),rpar(:,81))*ymax(ETAR) - y(ETAR))/tau(ETAR); 
dydt(DAG) = (OR(act(y(ETAR),rpar(:,152)),act(y(AT1R),rpar(:,153)))*ymax(DAG) - y(DAG))/tau(DAG); 
dydt(PKC) = (act(y(DAG),rpar(:,157))*ymax(PKC) - y(PKC))/tau(PKC); 
dydt(TRPC) = (act(y(DAG),rpar(:,154))*ymax(TRPC) - y(TRPC))/tau(TRPC); 
dydt(NE) = (act(y(inputNE),rpar(:,18))*ymax(NE) - y(NE))/tau(NE); 
dydt(inputNE) = (rpar(1,7)*ymax(inputNE) - y(inputNE))/tau(inputNE); 
dydt(BAR) = (act(y(NE),rpar(:,80))*ymax(BAR) - y(BAR))/tau(BAR); 
dydt(Forskolin) = (act(y(inputForskolin),rpar(:,22))*ymax(Forskolin) - y(Forskolin))/tau(Forskolin); 
dydt(inputForskolin) = (rpar(1,11)*ymax(inputForskolin) - y(inputForskolin))/tau(inputForskolin); 
dydt(AC) = (OR(act(y(BAR),rpar(:,93)),OR(AND(rpar(:,94),act(y(AT1R),rpar(:,94)),act(y(BAR),rpar(:,94))),act(y(Forskolin),rpar(:,137))))*ymax(AC) - y(AC))/tau(AC); 
dydt(cAMP) = (act(y(AC),rpar(:,95))*ymax(cAMP) - y(cAMP))/tau(cAMP); 
dydt(PKA) = (act(y(cAMP),rpar(:,63))*ymax(PKA) - y(PKA))/tau(PKA); 
dydt(CREB) = (act(y(PKA),rpar(:,78))*ymax(CREB) - y(CREB))/tau(CREB); 
dydt(CBP) = (OR(inhib(y(smad3),rpar(:,65)),inhib(y(CREB),rpar(:,66)))*ymax(CBP) - y(CBP))/tau(CBP); 
dydt(TGFB) = (OR(act(y(inputTGFB),rpar(:,13)),OR(AND(rpar(:,23),act(y(TGFBmRNA),rpar(:,23)),act(y(MMP9),rpar(:,23))),OR(AND(rpar(:,24),act(y(TGFBmRNA),rpar(:,24)),act(y(MMP2),rpar(:,24))),act(y(TGFBmRNA),rpar(:,25)))))*ymax(TGFB) - y(TGFB))/tau(TGFB); 
dydt(inputTGFB) = (rpar(1,2)*ymax(inputTGFB) - y(inputTGFB))/tau(inputTGFB); 
dydt(TGFB1R) = (AND(rpar(:,76),act(y(TGFB),rpar(:,76)),inhib(y(BAMBI),rpar(:,76)))*ymax(TGFB1R) - y(TGFB1R))/tau(TGFB1R); 
dydt(smad3) = (AND(rpar(:,47),act(y(TGFB1R),rpar(:,47)),inhib(y(smad7),rpar(:,47)),inhib(y(PKG),rpar(:,47)))*ymax(smad3) - y(smad3))/tau(smad3); 
dydt(smad7) = (act(y(STAT),rpar(:,138))*ymax(smad7) - y(smad7))/tau(smad7); 
dydt(TGFBmRNA) = (act(y(AP1),rpar(:,97))*ymax(TGFBmRNA) - y(TGFBmRNA))/tau(TGFBmRNA); 
dydt(BAMBI) = (AND(rpar(:,136),act(y(TGFB),rpar(:,136)),act(y(IL1RI),rpar(:,136)))*ymax(BAMBI) - y(BAMBI))/tau(BAMBI); 
dydt(PDGF) = (act(y(inputPDGF),rpar(:,19))*ymax(PDGF) - y(PDGF))/tau(PDGF); 
dydt(inputPDGF) = (rpar(1,8)*ymax(inputPDGF) - y(inputPDGF))/tau(inputPDGF); 
dydt(PDGFR) = (act(y(PDGF),rpar(:,92))*ymax(PDGFR) - y(PDGFR))/tau(PDGFR); 
dydt(NP) = (act(y(inputNP),rpar(:,21))*ymax(NP) - y(NP))/tau(NP); 
dydt(inputNP) = (rpar(1,10)*ymax(inputNP) - y(inputNP))/tau(inputNP); 
dydt(NPRA) = (act(y(NP),rpar(:,103))*ymax(NPRA) - y(NPRA))/tau(NPRA); 
dydt(cGMP) = (act(y(NPRA),rpar(:,104))*ymax(cGMP) - y(cGMP))/tau(cGMP); 
dydt(PKG) = (act(y(cGMP),rpar(:,105))*ymax(PKG) - y(PKG))/tau(PKG); 
dydt(mechanical) = (act(y(inputMechanical),rpar(:,14))*ymax(mechanical) - y(mechanical))/tau(mechanical); 
dydt(inputMechanical) = (rpar(1,3)*ymax(inputMechanical) - y(inputMechanical))/tau(inputMechanical); 
dydt(B1int) = (OR(AND(rpar(:,62),act(y(PKC),rpar(:,62)),act(y(mechanical),rpar(:,62))),act(y(mechanical),rpar(:,67)))*ymax(B1int) - y(B1int))/tau(B1int); 
dydt(Rho) = (act(y(B1int),rpar(:,43))*ymax(Rho) - y(Rho))/tau(Rho); 
dydt(ROCK) = (act(y(Rho),rpar(:,99))*ymax(ROCK) - y(ROCK))/tau(ROCK); 
dydt(Ca) = (act(y(TRPC),rpar(:,155))*ymax(Ca) - y(Ca))/tau(Ca); 
dydt(calcineurin) = (act(y(Ca),rpar(:,156))*ymax(calcineurin) - y(calcineurin))/tau(calcineurin); 
dydt(NFAT) = (act(y(calcineurin),rpar(:,147))*ymax(NFAT) - y(NFAT))/tau(NFAT); 
dydt(IL6) = (OR(act(y(inputIL6),rpar(:,15)),act(y(IL6mRNA),rpar(:,31)))*ymax(IL6) - y(IL6))/tau(IL6); 
dydt(inputIL6) = (rpar(1,4)*ymax(inputIL6) - y(inputIL6))/tau(inputIL6); 
dydt(IL6mRNA) = (OR(AND(rpar(:,28),act(y(CREB),rpar(:,28)),act(y(CBP),rpar(:,28))),OR(act(y(NFKB),rpar(:,29)),act(y(AP1),rpar(:,30))))*ymax(IL6mRNA) - y(IL6mRNA))/tau(IL6mRNA); 
dydt(gp130) = (act(y(IL6),rpar(:,37))*ymax(gp130) - y(gp130))/tau(gp130); 
dydt(STAT) = (act(y(gp130),rpar(:,44))*ymax(STAT) - y(STAT))/tau(STAT); 
dydt(IL1) = (act(y(inputIL1),rpar(:,16))*ymax(IL1) - y(IL1))/tau(IL1); 
dydt(inputIL1) = (rpar(1,5)*ymax(inputIL1) - y(inputIL1))/tau(inputIL1); 
dydt(IL1RI) = (act(y(IL1),rpar(:,83))*ymax(IL1RI) - y(IL1RI))/tau(IL1RI); 
dydt(TNFa) = (act(y(inputTNFa),rpar(:,17))*ymax(TNFa) - y(TNFa))/tau(TNFa); 
dydt(inputTNFa) = (rpar(1,6)*ymax(inputTNFa) - y(inputTNFa))/tau(inputTNFa); 
dydt(TNFaR) = (act(y(TNFa),rpar(:,102))*ymax(TNFaR) - y(TNFaR))/tau(TNFaR); 
dydt(NFKB) = (OR(act(y(IL1RI),rpar(:,41)),OR(act(y(ERK),rpar(:,53)),OR(act(y(p38),rpar(:,54)),act(y(Akt),rpar(:,133)))))*ymax(NFKB) - y(NFKB))/tau(NFKB); 
dydt(PI3K) = (OR(act(y(TNFaR),rpar(:,45)),OR(act(y(TGFB1R),rpar(:,129)),OR(act(y(PDGFR),rpar(:,130)),act(y(FAK),rpar(:,131)))))*ymax(PI3K) - y(PI3K))/tau(PI3K); 
dydt(Akt) = (act(y(PI3K),rpar(:,128))*ymax(Akt) - y(Akt))/tau(Akt); 
dydt(p38) = (OR(act(y(ROS),rpar(:,39)),OR(act(y(MKK3),rpar(:,110)),OR(act(y(Ras),rpar(:,127)),AND(rpar(:,140),act(y(Rho),rpar(:,140)),inhib(y(Rac1),rpar(:,140))))))*ymax(p38) - y(p38))/tau(p38); 
dydt(TRAF) = (OR(act(y(TGFB1R),rpar(:,111)),act(y(TNFaR),rpar(:,119)))*ymax(TRAF) - y(TRAF))/tau(TRAF); 
dydt(ASK1) = (OR(act(y(TRAF),rpar(:,120)),act(y(IL1RI),rpar(:,123)))*ymax(ASK1) - y(ASK1))/tau(ASK1); 
dydt(MKK3) = (act(y(ASK1),rpar(:,121))*ymax(MKK3) - y(MKK3))/tau(MKK3); 
dydt(PP1) = (act(y(p38),rpar(:,109))*ymax(PP1) - y(PP1))/tau(PP1); 
dydt(JNK) = (OR(act(y(ROS),rpar(:,40)),OR(AND(rpar(:,114),inhib(y(NFKB),rpar(:,114)),act(y(MKK4),rpar(:,114))),AND(rpar(:,141),inhib(y(Rho),rpar(:,141)),act(y(MKK4),rpar(:,141)))))*ymax(JNK) - y(JNK))/tau(JNK); 
dydt(abl) = (act(y(PDGFR),rpar(:,115))*ymax(abl) - y(abl))/tau(abl); 
dydt(Rac1) = (OR(act(y(B1int),rpar(:,42)),act(y(abl),rpar(:,116)))*ymax(Rac1) - y(Rac1))/tau(Rac1); 
dydt(MEKK1) = (OR(act(y(FAK),rpar(:,96)),act(y(Rac1),rpar(:,112)))*ymax(MEKK1) - y(MEKK1))/tau(MEKK1); 
dydt(MKK4) = (OR(act(y(MEKK1),rpar(:,113)),act(y(ASK1),rpar(:,122)))*ymax(MKK4) - y(MKK4))/tau(MKK4); 
dydt(ERK) = (act(y(MEK1),rpar(:,108))*ymax(ERK) - y(ERK))/tau(ERK); 
dydt(Ras) = (act(y(AT1R),rpar(:,148))*ymax(Ras) - y(Ras))/tau(Ras); 
dydt(Raf) = (act(y(Ras),rpar(:,106))*ymax(Raf) - y(Raf))/tau(Raf); 
dydt(MEK1) = (OR(AND(rpar(:,38),act(y(ROS),rpar(:,38)),inhib(y(PP1),rpar(:,38))),AND(rpar(:,107),inhib(y(ERK),rpar(:,107)),act(y(Raf),rpar(:,107))))*ymax(MEK1) - y(MEK1))/tau(MEK1); 
dydt(FAK) = (act(y(ROCK),rpar(:,149))*ymax(FAK) - y(FAK))/tau(FAK); 
dydt(epac) = (act(y(cAMP),rpar(:,98))*ymax(epac) - y(epac))/tau(epac); 
dydt(Factin) = (act(y(ROCK),rpar(:,142))*ymax(Factin) - y(Factin))/tau(Factin); 
dydt(FA) = (AND(rpar(:,144),act(y(B1int),rpar(:,144)),act(y(Factin),rpar(:,144)))*ymax(FA) - y(FA))/tau(FA); 
dydt(migration) = (OR(act(y(MMP9),rpar(:,100)),OR(act(y(MMP2),rpar(:,101)),AND(rpar(:,125),inhib(y(PKA),rpar(:,125)),act(y(epac),rpar(:,125)))))*ymax(migration) - y(migration))/tau(migration); 
dydt(cmyc) = (act(y(JNK),rpar(:,117))*ymax(cmyc) - y(cmyc))/tau(cmyc); 
dydt(CTGF) = (AND(rpar(:,48),act(y(CBP),rpar(:,48)),act(y(smad3),rpar(:,48)),act(y(ERK),rpar(:,48)))*ymax(CTGF) - y(CTGF))/tau(CTGF); 
dydt(proliferation) = (OR(act(y(AP1),rpar(:,77)),OR(act(y(CREB),rpar(:,79)),OR(act(y(CTGF),rpar(:,82)),OR(act(y(PKC),rpar(:,84)),act(y(cmyc),rpar(:,118))))))*ymax(proliferation) - y(proliferation))/tau(proliferation); 
dydt(SRF) = (act(y(Factin),rpar(:,143))*ymax(SRF) - y(SRF))/tau(SRF); 
dydt(EDAFN) = (OR(AND(rpar(:,64),act(y(CBP),rpar(:,64)),act(y(smad3),rpar(:,64))),OR(act(y(NFAT),rpar(:,68)),act(y(NFKB),rpar(:,134))))*ymax(EDAFN) - y(EDAFN))/tau(EDAFN); 
dydt(aSMA) = (OR(AND(rpar(:,146),act(y(CBP),rpar(:,146)),act(y(smad3),rpar(:,146)),act(y(SRF),rpar(:,146))),OR(AND(rpar(:,150),act(y(CBP),rpar(:,150)),act(y(smad3),rpar(:,150))),act(y(SRF),rpar(:,151))))*ymax(aSMA) - y(aSMA))/tau(aSMA); 
dydt(AP1) = (OR(act(y(ERK),rpar(:,57)),act(y(JNK),rpar(:,135)))*ymax(AP1) - y(AP1))/tau(AP1); 
dydt(TIMP1) = (act(y(AP1),rpar(:,60))*ymax(TIMP1) - y(TIMP1))/tau(TIMP1); 
dydt(TIMP2) = (act(y(AP1),rpar(:,61))*ymax(TIMP2) - y(TIMP2))/tau(TIMP2); 
dydt(PAI1) = (act(y(smad3),rpar(:,124))*ymax(PAI1) - y(PAI1))/tau(PAI1); 
dydt(proMMP14) = (OR(act(y(AP1),rpar(:,91)),act(y(NFKB),rpar(:,126)))*ymax(proMMP14) - y(proMMP14))/tau(proMMP14); 
dydt(proMMP1) = (OR(AND(rpar(:,55),inhib(y(smad3),rpar(:,55)),act(y(NFKB),rpar(:,55))),AND(rpar(:,160),inhib(y(smad3),rpar(:,160)),act(y(AP1),rpar(:,160))))*ymax(proMMP1) - y(proMMP1))/tau(proMMP1); 
dydt(proMMP2) = (OR(act(y(STAT),rpar(:,49)),act(y(AP1),rpar(:,58)))*ymax(proMMP2) - y(proMMP2))/tau(proMMP2); 
dydt(proMMP9) = (OR(act(y(STAT),rpar(:,50)),AND(rpar(:,59),act(y(NFKB),rpar(:,59)),act(y(AP1),rpar(:,59))))*ymax(proMMP9) - y(proMMP9))/tau(proMMP9); 
dydt(MMP1) = (AND(rpar(:,72),inhib(y(TIMP1),rpar(:,72)),act(y(proMMP1),rpar(:,72)))*ymax(MMP1) - y(MMP1))/tau(MMP1); 
dydt(MMP2) = (OR(AND(rpar(:,74),inhib(y(TIMP1),rpar(:,74)),act(y(proMMP2),rpar(:,74)),act(y(MMP14),rpar(:,74))),AND(rpar(:,75),act(y(proMMP2),rpar(:,75)),act(y(MMP14),rpar(:,75))))*ymax(MMP2) - y(MMP2))/tau(MMP2); 
dydt(MMP9) = (OR(AND(rpar(:,71),inhib(y(TIMP1),rpar(:,71)),act(y(proMMP9),rpar(:,71))),AND(rpar(:,73),inhib(y(TIMP2),rpar(:,73)),act(y(proMMP9),rpar(:,73))))*ymax(MMP9) - y(MMP9))/tau(MMP9); 
dydt(MMP14) = (act(y(proMMP14),rpar(:,70))*ymax(MMP14) - y(MMP14))/tau(MMP14); 
dydt(periostin) = (OR(AND(rpar(:,51),act(y(CBP),rpar(:,51)),act(y(smad3),rpar(:,51))),AND(rpar(:,52),act(y(CREB),rpar(:,52)),act(y(CBP),rpar(:,52))))*ymax(periostin) - y(periostin))/tau(periostin); 
dydt(CImRNA) = (OR(AND(rpar(:,85),act(y(CBP),rpar(:,85)),act(y(smad3),rpar(:,85)),inhib(y(epac),rpar(:,85))),act(y(SRF),rpar(:,139)))*ymax(CImRNA) - y(CImRNA))/tau(CImRNA); 
dydt(CIIImRNA) = (OR(AND(rpar(:,86),act(y(CBP),rpar(:,86)),act(y(smad3),rpar(:,86)),inhib(y(epac),rpar(:,86))),act(y(SRF),rpar(:,145)))*ymax(CIIImRNA) - y(CIIImRNA))/tau(CIIImRNA); 
dydt(CI) = (OR(AND(rpar(:,87),inhib(y(MMP1),rpar(:,87)),act(y(CImRNA),rpar(:,87))),AND(rpar(:,89),inhib(y(MMP2),rpar(:,89)),act(y(CImRNA),rpar(:,89))))*ymax(CI) - y(CI))/tau(CI); 
dydt(CIII) = (OR(AND(rpar(:,88),inhib(y(MMP1),rpar(:,88)),act(y(CIIImRNA),rpar(:,88))),AND(rpar(:,90),inhib(y(MMP2),rpar(:,90)),act(y(CIIImRNA),rpar(:,90))))*ymax(CIII) - y(CIII))/tau(CIII); 
dydt(LOX) = (act(y(Akt),rpar(:,158))*ymax(LOX) - y(LOX))/tau(LOX); 
dydt(xlinkFibers) = (AND(rpar(:,159),act(y(CI),rpar(:,159)),act(y(LOX),rpar(:,159)))*ymax(xlinkFibers) - y(xlinkFibers))/tau(xlinkFibers); 
dydt(contraction) = (OR(act(y(FA),rpar(:,161)),act(y(aSMA),rpar(:,162)))*ymax(contraction) - y(contraction))/tau(contraction); 

dydt = dydt.*KI';
end

% utility functions 
function calcFact = calculateFact(EC50, n, w, x)
    beta = (EC50.^n - 1)./(2*EC50.^n - 1); 
    K = (beta - 1).^(1./n); 
    calcFact = w.*(beta.*x.^n)./(K.^n + x.^n);
    
    if calcFact > w 
        calcFact = w; % cap fact(x)<= 1 
    elseif calcFact > 1
        calcFact = 1;
    elseif calcFact < 0
        calcFact = 0;
    end  
end

function fact = act(x,rpar) 
% hill activation function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    dose = rpar(4);
    drugBinding = rpar(5);
    drugAgonism = rpar(6);
    
    if drugBinding == 1 % Competitive drug
        x = x - dose;
    end
    if x < 0
        x = 0;
    end
    if x > 1
        x = 1;
    end
    
    fact = calculateFact(EC50, n, w, x);
    if drugBinding == -1 && drugAgonism == 1 % Non-Competitive Agonist
        fact = fact * (1-dose) + dose;
    elseif drugBinding == -1 && drugAgonism == -1 % Non-Competitive Antagonist
        fact = fact * (1-dose);
    end
end

function finhib = inhib(x,rpar) 
% inverse hill function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    dose = rpar(4);
    drugBinding = rpar(5);
    drugAgonism = rpar(6);
    
    if drugBinding == 1 % Competitive drug
        x = x + dose; % THIS IS THE OPPOSITE ('+', '-') FOR INHIBITION RXNS
    end
    if x < 0
        x = 0;
    end
    if x > 1
        x = 1;
    end
    
    fact = calculateFact(EC50, n, w, x);
    finhib = w - fact;
    
    if drugBinding == -1 && drugAgonism == 1 % Non-Competitive Agonist
        finhib = finhib * (1-dose) + dose;
    elseif drugBinding == -1 && drugAgonism == -1 % Non-Competitive Antagonist
        finhib = finhib * (1-dose);
    end
end

function z = OR(x,y) 
% OR logic gate 
    z = x + y - x*y;
end

function z = AND(rpar,varargin) 
% AND logic gate, multiplying all of the reactants together 
    w = rpar(1); 

    if w == 0 
        z = 0; 
    else 
        v = cell2mat(varargin); 
        z = prod(v)/w^(nargin-2);  
    end 
end