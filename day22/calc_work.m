clear all
c=constants;
wtA=14.e-3;
pressA=900.e2;
tempA=25 + c.Tc;

TdA=findTdwv(wtA,pressA);
thetaeA=thetaep(TdA,tempA,pressA);
wtB=wtA;
pressB=700.e2;
TdB=findTdwv(wtB,pressB);
thetaeB=thetaeA;
[tempB,wvB,wlB]=tinvert_thetae(thetaeB, wtB, pressB);

wtC=wtA;
pressC=900.e2;
TdC=findTdwv(wtC,pressC);
tempC=tempB;
thetaeC=thetaep(TdC,tempC,pressC);
skew=30.;
figHandle=figure(1);
[figureHandle,outputws,handlews]=makeSkew(figHandle,skew);
xtempA=convertTempToSkew(tempA - c.Tc,pressA*0.01,skew);
xtempB=convertTempToSkew(tempB - c.Tc,pressB*0.01,skew);
xtempC=convertTempToSkew(tempC - c.Tc,pressC*0.01,skew);
text(xtempA,pressA*0.01,'A',...,
            'edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempB,pressB*0.01,'B',...
            'edgecolor','b','fontweight','bold','fontsize',22,'color','b');
text(xtempC,pressC*0.01,'C',...
            'edgecolor','b','fontweight','bold','fontsize',22,'color','b');
pressLevs=linspace(700,900,60)*100.;
for i=1:numel(pressLevs)
    thePress=pressLevs(i);
    [temp,wv,wl]=tinvert_thetae(thetaeA, wtA, thePress);
    lineAB(i)=temp;
    rho=thePress/(c.Rd*temp);
    rhoAB(i)=rho;
end
tempCA=linspace(tempC,tempA,100);
rhoCA=pressA./(c.Rd*tempCA);
press900Vec=NaN(size(rhoCA));
press900Vec(:)=pressA;
rhoBC=pressLevs/(c.Rd*tempB);
xtemp=convertTempToSkew(lineAB - c.Tc,pressLevs*0.01,skew);    
semilogy(xtemp,pressLevs*0.01,'k-.','linewidth',2);
semilogy([xtempB,xtempC],[700.,900.],'b-.','linewidth',2);
semilogy([xtempC,xtempA],[900.,900.],'r-.','linewidth',2);
title('heat engine problem');
xleft=convertTempToSkew(10,1000.,skew);
xright=convertTempToSkew(30,1000.,skew);
axis([xleft,xright,650,1000.]);
print -depsc prob1.eps

figure(2)
clf;
plot(1./rhoCA,press900Vec*1.e-2);
ylim([700,1000.]);
hold on;
plot(1./rhoAB,pressLevs*1.e-2);
plot(1./rhoBC,pressLevs*1.e-2);
title('volume - pressure plot')
hold off;

alphaAB = 1./rhoAB;
alphaCA = 1./rhoCA;
alphaBC = 1./rhoBC;

%work done (in J) by air from A -> B  (should be positive, air is expanding)
w_AB = quad(@(alpha_interp) interp1(alphaAB, pressLevs, alpha_interp), alphaAB(end), alphaAB(1))

%work done by air from B -> C (should be negative, air is being compressed)
w_BC = quad(@(alpha_interp) interp1(alphaBC, pressLevs, alpha_interp), alphaBC(1) , alphaBC(end))

%work done by air from C -> A (isobaric process)
%(should be positive; work is being done on surroundings since there is heat flow 
% into the system and the pressure is constant)
w_CA = quad(@(alpha_interp) interp1(alphaCA, press900Vec, alpha_interp), alphaCA(1) , alphaCA(end))
%w_CA = p_A.*(alphaA - alphaC)

%total work done by air on surroundings (J)
totalwork = w_AB + w_BC + w_CA
