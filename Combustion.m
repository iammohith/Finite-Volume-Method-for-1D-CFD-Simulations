s=0.1; %stroke in meters
r=s/2; % radius in meters
l=4*s; % connecting rod length in metrs
tmax=3000; %peak temperature in K
rc=9; %compression ratio
d=s; %diameter in meter
cv=0.71; %the molar specific heat at constant volume
gamma=1.4; %the heat capacity ratio
n=1000; %r.p.m
theta_deg=-180:1:180; %cranck angle
theta=pi/180*theta_deg; %cranck angle in radians
p=101325; %inlet pressure
t=25; %inlet temperature in degC
t=t+273; %converting degC to K
R=0.287; %air gas constant
%% space conservation for constant gamma
y=zeros(1,361); %pre allocating array for position
ys=zeros(1,361); %pre allocating array for velocity
ya=zeros(1,361); %pre allocating array for acceleration
vs=zeros(1,362); %pre allocating array for swept volume
ts=zeros(1,362); %pre allocating array for instant temperature
ps=zeros(1,362); %pre allocating array for instant pressure
wcom=zeros(1,361); %pre allocating array for compression work
wexp=zeros(1,361); %pre allocating array for expansion work
%% space conservation for variable gamma
gammac=zeros(1,362); %pre allocating array for gamma as a function in temperature
cvc=zeros(1,362); %pre allocating array for cv as a function in temperature
tsc=zeros(1,362); %pre allocating array for instant temperature due to change of gamma
psc=zeros(1,362); %pre allocating array for instant pressure due to change of gamma
g=zeros(1,20); %pre allocating array for variable gamma
error=zeros(1,361); %pre allocating array for error in gamma
%% Calculated constants
w=2*pi*n/60;
y(1)=l*((r/l)^2*sin(theta(1)/2)^2)+r*(1-cos(theta(1)));
ys(1)=r*w*(sin(theta(1))+r/l*sin(2*theta(1))/2);
ya(1)=r*w^2*(cos(theta(1))+r/l*cos(2*theta(1)));
vmax=pi/4*d^2*y(1);% max swept volume
vc=vmax/(rc-1); %clearance volume
vs(1)=vmax+vc;
ts(1)=t;
ps(1)=p;
tsc(1)=t;
psc(1)=p;
gammac(1)=1.4;
cvc(1)=0.71;
vs(362)=vmax+vc;
ts(362)=t;
ps(362)=p;
tsc(362)=t;
psc(362)=p;
%% Equations and loops
for i=2:361
    if i<181
        
        y(i)=l*((r/l)^2*sin(theta(i)/2)^2)+r*(1-cos(theta(i)));
        ys(i)=r*w*(sin(theta(i))+r/l*sin(2*theta(i))/2);
        ya(i)=r*w^2*(cos(theta(i))+r/l*cos(2*theta(i)));
        vs(i)=vc+pi/4*d^2*y(i);
        
        %% for contant gamma
        ps(i)=(vs(i-1)/vs(i))^gamma*ps(i-1);
        ts(i)=(vs(i-1)/vs(i))^(gamma-1)*ts(i-1);
        
        
        %% for variable gamma
        g(1)=gammac(i-1);
        
        for j=2:20
            t2=(vs(i-1)/vs(i))^(g(j-1)-1)*tsc(i-1);
            cv2=(0.71*(t2-tsc(i-1))+19*10^(-5)/2*((t2)^2-(tsc(i-1))^2))/(t2-tsc(i-1));
            g(j)=(cv2+R)/cv2;
            error(i)=(abs((g(j-1)-g(j))))/g(j);
        end
        
        gammac(i)=g(20);
        cvc(i)=cv2;
        tsc(i)=t2;
        psc(i)=(vs(i-1)/vs(i))^gammac(i)*psc(i-1);
        
    elseif i==181
        
        y(i)=l*((r/l)^2*sin(theta(i)/2)^2)+r*(1-cos(theta(i)));
        ys(i)=r*w*(sin(theta(i))+r/l*sin(2*theta(i))/2);
        ya(i)=r*w^2*(cos(theta(i))+r/l*cos(2*theta(i)));
        vs(i)=vc;
        %% for constant gamma
        ts(i)=tmax;
        ps(i)=(ts(i)/ts(i-1))*ps(i-1);
        
        %% for variable gamma
        tsc(i)=tmax;
        psc(i)=(tsc(i)/tsc(i-1))*psc(i-1);
        cvc(i)=0.71+19*10^(-5)/2*(tsc(i)-tsc(i-1));
        gammac(i)=(cvc(i)+R)/cvc(i);
        
        
    else
        y(i)=l*((r/l)^2*sin(theta(i)/2)^2)+r*(1-cos(theta(i)));
        ys(i)=r*w*(sin(theta(i))+r/l*sin(2*theta(i))/2);
        ya(i)=r*w^2*(cos(theta(i))+r/l*cos(2*theta(i)));
        vs(i)=vc+pi/4*d^2*y(i);
        
        %% for constant gamma
        ps(i)=(vs(i-1)/vs(i))^gamma*ps(i-1);
        ts(i)=(vs(i-1)/vs(i))^(gamma-1)*ts(i-1);
        
        
        %% for variable gamma
        g(1)=gammac(i-1);
        
        for j=2:20
            t2=(vs(i-1)/vs(i))^(g(j-1)-1)*tsc(i-1);
            cv2=(0.71*(t2-tsc(i-1))+19*10^(-5)/2*((t2)^2-(tsc(i-1))^2))/(t2-tsc(i-1));
            g(j)=(cv2+R)/cv2;
            error(i)=(abs((g(j-1)-g(j))))/g(j);
        end
        
        gammac(i)=g(20);
        cvc(i)=cv2;
        tsc(i)=t2;
        psc(i)=(vs(i-1)/vs(i))^gammac(i)*psc(i-1);
        
    end
    
end
%% from data calculated above for constant gamma
qrej=cv*(ts(361)-ts(362));% q rejected
qadded=cv*(ts(181)-ts(180));% q added
wcomt=cv*(ts(180)-ts(1));%total compression work
wexpt=cv*(ts(181)-ts(361));%total expansion work
wt=wexpt-wcomt;%total cycle work
eta=wt/qadded; %efficiency
%% from data calculated above for variable gamma
gammac(362)=1.4;
cvc(362)=0.71+19*10^(-5)/2*(tsc(361)+tsc(362));
qaddedc=(0.71+19*10^(-5)/2*(tsc(180)+tsc(181)))*(tsc(181)-tsc(180));
qrejc=cvc(362)*(tsc(361)-tsc(362));
wcomtc=(0.71+19*10^(-5)/2*(tsc(1)+tsc(180)))*(tsc(180)-tsc(1));
wexptc=(0.71+19*10^(-5)/2*(tsc(180)+tsc(361)))*(tsc(181)-tsc(361));
wtc=wexptc-wcomtc;
etac=wtc/qaddedc;
%% figures
figure;
plot(theta_deg,y,'Color',[0.8 0.3 0.01],'LineWidth',1.5,'LineStyle','-','Marker','none','MarkerSize',15)
title('change of displacement with crank angle')
xlabel('crank angle in degrees')
ylabel('piston displacement(m)')
ax = gca; ax.FontName = 'Century Gothic'; ax.FontSize=16; ax.FontWeight='bold'; ax.GridAlpha=0.07; ax.GridLineStyle='--'; ax.LineWidth=1.5; ax.XColor=[0 0 0]; ax.XMinorTick='on'; ax.YColor=[0 0 0]; ax.YMinorTick='on';
figure;
plot(theta_deg,ys,'Color',[0.8 0.3 0.01],'LineWidth',1.5,'LineStyle','-','Marker','none','MarkerSize',15)
xlabel('crank angle in degrees')
ylabel('piston velocity(m/sec)')
title('change of piston velocity with crank angle')
ax = gca; ax.FontName = 'Century Gothic'; ax.FontSize=16; ax.FontWeight='bold'; ax.GridAlpha=0.07; ax.GridLineStyle='--'; ax.LineWidth=1.5; ax.XColor=[0 0 0]; ax.XMinorTick='on'; ax.YColor=[0 0 0]; ax.YMinorTick='on';
figure;
plot(theta_deg,ya,'Color',[0.8 0.3 0.01],'LineWidth',1.5,'LineStyle','-','Marker','none','MarkerSize',15)
xlabel('crank angle in degrees')
ylabel('piston accceleration(m/sec^2)')
title('change of piston displacement with crank angle')
ax = gca; ax.FontName = 'Century Gothic'; ax.FontSize=16; ax.FontWeight='bold'; ax.GridAlpha=0.07; ax.GridLineStyle='--'; ax.LineWidth=1.5; ax.XColor=[0 0 0]; ax.XMinorTick='on'; ax.YColor=[0 0 0]; ax.YMinorTick='on';
figure;
plot(vs,ps,'b',vs,psc,'r','LineWidth',1.5,'LineStyle','-','Marker','none','MarkerSize',15)
legend('temp independent properties ','temp dependent properties')
xlabel('volume (m^3/sec)')
ylabel('pressure(Pa)')
title('cycle PV diagram')
ax = gca; ax.FontName = 'Century Gothic'; ax.FontSize=16; ax.FontWeight='bold'; ax.GridAlpha=0.07; ax.GridLineStyle='--'; ax.LineWidth=1.5; ax.XColor=[0 0 0]; ax.XMinorTick='on'; ax.YColor=[0 0 0]; ax.YMinorTick='on';
figure;
plot(theta_deg,ps(1:361),'b',theta_deg,psc(1:361),'r','LineWidth',1.5,'LineStyle','-','Marker','none','MarkerSize',15)
legend('temp independent properties ','temp dependent properties')
xlabel('crank angle in degrees')
ylabel('pressure(Pa)')
title('change of pressure with crank angle')
ax = gca; ax.FontName = 'Century Gothic'; ax.FontSize=16; ax.FontWeight='bold'; ax.GridAlpha=0.07; ax.GridLineStyle='--'; ax.LineWidth=1.5; ax.XColor=[0 0 0]; ax.XMinorTick='on'; ax.YColor=[0 0 0]; ax.YMinorTick='on';
figure;
plot(gammac(1:180),tsc(1:180),'c',gammac(180:181),tsc(180:181),'r',gammac(181:361),tsc(181:361),'b','LineWidth',1.5,'LineStyle','-','Marker','none','MarkerSize',15)
legend('compression','combustion','expansion')
xlabel('gamma')
ylabel('temperature(k)')
title('change of gamma with temperature')
ax = gca; ax.FontName = 'Century Gothic'; ax.FontSize=16; ax.FontWeight='bold'; ax.GridAlpha=0.07; ax.GridLineStyle='--'; ax.LineWidth=1.5; ax.XColor=[0 0 0]; ax.XMinorTick='on'; ax.YColor=[0 0 0]; ax.YMinorTick='on';