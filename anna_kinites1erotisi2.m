clc
clear all
clf
h=20; %ypsos ebodiou
ht=30; %ypsos pompou
hr=1.8; %ypsos dekti
d1=500; %apostasi metaxi empodiou kai pompou
d2= 510:5:1400; %apostasi metaxi empodiou kai dekti
epsilon1 = 8.854187*(10.^(-12));
sigma=0.01; 
EpsilonRC=epsilon1-j*sigma;
lamda=0.333;
R1=sqrt((lamda*(d1*d2))/(d1+d2));
%%%%%%%%%%%%%%%%%%%%%%%%% Geometrikoi parametroi %%%%%%%%%%%%%%%%%%%%%
TR = (sqrt(d1+d2).^2+(ht-hr)^2);
h1 = h-((ht*d2)+(hr*d1)./(d1+d2));
TR2 = sqrt((d1+d2).^2+(ht+hr).^2); %T'R
h2 = h+(((h*d2)-(hr*d1))./(d1+d2));
h3 = h-(((ht*d2)-(hr*d1))./d1+d2);
h4 = h+((ht*d2)+(hr*d1)./(d1+d2));

%%%%%%%%%%%%%%%%%%%%%%%% syntelestes anaklasis %%%%%%%%%%%%%%%%%%%%%%%%
a= atand((ht+h)./d1);
b= atand((hr+h)./d2);

Rv=(sin(a)-sqrt(EpsilonRC-(cos(a).^2)))/(sin(a)+sqrt((EpsilonRC-(cos(a).^2))));  %R vertical
Rh=((EpsilonRC*sin(a))-sqrt(EpsilonRC-(cos(a).^2)))/((EpsilonRC*sin(a))+sqrt((EpsilonRC-(cos(a).^2)))); %R horizontal
Ga=(Rv+Rh)/(1+(Rv*Rh));

Rvb=(sin(b)-sqrt(EpsilonRC-(cos(b).^2)))/(sin(b)+sqrt((EpsilonRC-(cos(b).^2))));  %R vertical
Rhb=((EpsilonRC*sin(b))-sqrt(EpsilonRC-(cos(b).^2)))/((EpsilonRC*sin(b))+sqrt((EpsilonRC-(cos(b).^2)))); %R horizontal
Gb=(Rvb+Rhb)/(1+(Rvb*Rhb));

%%%%%%%%%%%%%%%%%%%%%%    free space field   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E01=exp(-j*2*pi/lamda*TR)/TR;
E02=exp(-j*2*pi/lamda*TR2)/TR2;
E03=E02;
E04=E01;

%%%%%%%%%%%%%%%%%%%%%%%% FRESNEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1=sqrt(2)*(h1/R1);
u2=sqrt(2)*(h2/R1);
u3=sqrt(2)*(h3/R1);
u4=sqrt(2)*(h4/R1);

%%%%%%%ypologismos u1%%%%%%%%%%%%%%%%
Cu1=mfun('FresnelC', u1);
Su1=mfun('FresnelS', u1);
F1=Cu1-j*Su1;
%%%%u2
Cu2=mfun('FresnelC', u2);
Su2=mfun('FresnelS', u2);
F2=Cu2-j*Su2;
%%%%%u3
Cu3=mfun('FresnelC', u3);
Su3=mfun('FresnelS', u3);
F3=Cu3-j*Su3;
%%%%%%u4
Cu4=mfun('FresnelC', u4);
Su4=mfun('FresnelS', u4);
F4=Cu4-j*Su4;



Er_E01 = 0.5*(1-(1+j).*F1); %Volt
Er_E02 = 0.5*(1-(1+j).*F2)*Ga;
Er_E03 = 0.5*(1-(1+j).*F3)*Gb;
Er_E04 = 0.5*(1-(1+j).*F4)*Ga*Gb;

Etotal=(Er_E01*E01)+(Er_E02*E02)+(Er_E03*E03)+(Er_E04*E04);   %%% to synoliko pedio ton 4 synistoswn kai tou apey8eias pediou

plot(d1+d2,(20*log10(abs(Etotal))),'gx')        %%%platos
grid on
xlabel('Fresnel Parameters')
ylabel('')
figure ()
plot(d1+d2,angle(20*log10(Etotal)))      %%%%fasi
figure ()
plot(d2,(20*log10(abs(E01+E02+E03+E04./Etotal))),'m')    %%%apoleies
hold on

%%%%%%%%%%%%%%%%%%%%%%%plot ton epimerous simeion%%%%%%%%%%%%%%%%%%%
% plot(d1+d2,(20*log10(abs(Er_E01))),'gx')
% grid on
% xlabel('Fresnel Parameters')
% ylabel('')
% figure ()
% plot(d1+d2,angle(20*log10(Er_E01)))
% figure ()
% plot(d1+d2,(20*log10(abs(E01/Etotal)),'m')
% hold on
% 
% plot(d1+d2,(20*log10(abs(Er_E02))),'gx')
% grid on
% xlabel('Fresnel Parameters')
% ylabel('')
% figure ()
% plot(d1+d2,angle(20*log10(Er_E02)))
% figure ()
% plot(d1+d2,(10*log10(abs(E02/Etotal)),'m')
% hold on
% 
% plot(d1+d2,(20*log10(abs(Er_E03))),'gx')
% grid on
% xlabel('Fresnel Parameters')
% ylabel('')
% figure ()
% plot(d1+d2,angle(20*log10(Er_E03)))
% figure ()
% plot(d1+d2,(10*log10(abs(E03/Etotal)),'m')
% hold on
% 
% 
% plot(d1+d2,(20*log10(abs(Er_E04))),'gx'))
% grid on
% xlabel('Fresnel Parameters')
% ylabel('')
% figure ()
% plot(d1+d2,angle(20*log10(Er_E04)))
% figure ()
% plot(d1+d2,(10*log10(abs(E02/Etotal)),'m')
% hold on
