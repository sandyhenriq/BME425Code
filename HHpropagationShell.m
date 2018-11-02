
%function HHpropagation
% set constants of interest 
R=8.314;
F=96487;
T=25+273.15;
RTF=1000*(R*T/F);  %mV

Cm=1.0; % microF/cm^2
 
%set concentrations
Ke=20;
Ki=326;
Nae=381;
Nai=50;
%set conductances mS/cm^2
gk=36.;
gna=120.;
gl=0.3;

%calculate Nernst Potentials
Ena=RTF*log(Nae/Nai);
Ek=RTF*log(Ke/Ki);
Vrest=-60;

% Variable cnt for printing frequency for animation
cnt=0;

%defining parameters for axon
Numnodes=150;
radius = 0.0015; %cm
dx=0.01;
x=[0:dx:(Numnodes-1)*dx]; %Space Vector
Ri=0.5; %Axial Resistivity-kohms-cm
coeff = radius/(2*Ri*dx^2);
%resting  potential and gates values
Vcalc=Vrest.*[ones(1,Numnodes)];
Iion=[zeros(1,Numnodes)];
smallVm=Vcalc-Vrest;

%note vector math - this is specific to the way I wrote my functions
ngate=alphan(smallVm)./(alphan(smallVm)+betan(smallVm));
hgate=alphah(smallVm)./(alphah(smallVm)+betah(smallVm));
mgate=alpham(smallVm)./(alpham(smallVm)+betam(smallVm));

 %vector math
El = (gk*(ngate(1))*(Vrest-Ek)+gna*hgate(1)*(mgate(1)^3)*(Vrest-Ena)+gl*Vrest)/gl;

%stimulus settings
  Stimulus=100;
  tstart=0.1;
  tdur=0.5;
  
 
  dt=0.001;
  numsteps=20/dt;
  tvect=dt*[0:numsteps-1];  %Time Vector

  % printing Vm vs space at 3 times
  Vprint1=Vrest.*[ones(1,Numnodes)];
  Vprint2=Vrest.*[ones(1,Numnodes)];
  Vprint3=Vrest.*[ones(1,Numnodes)];
  time1 = 2;
  time2 = 4;
  time3 = 6;
  pcount=0;
  
  %Connectivity Matrix A
  %THIS IS THE PART YOU NEED TO FILL IN
  i = 1;
  j = 1;
  m = Numnodes;
  A = [];
  for i = 1:m
      for j = 1:m
          if i == 1 && j ==1
              A(i,j) = -1;
          elseif i == j
              A(i,j) = -2;
          elseif i==m && j==m
              A(i,j) = -1;
          elseif i-j == 1 || i-j == -1
              A(i,j) = 1;
          else
              A(i,j) = 0;
          end
      end
      A(m,m) = -1;
  end
  
  Istim=[zeros(1,Numnodes)];
  
  %beginning storing arrays for 3 nodes for time plots
  outvm1=[zeros(1,Numnodes)];
  outvm1(1)=Vcalc(1);
  outvm2=[zeros(1,Numnodes)];
  outvm2(1)=Vcalc(50);
  outvm3=[zeros(1,Numnodes)];
  outvm3(1)=Vcalc(100);
  
 for i=2:numsteps
     Istim=[zeros(1,Numnodes)];
     %apply stimulus to first 3 nodes if it is on
     if tvect(i)>tstart && tvect(i)<tstart+tdur
         Istim(1:3)=Stimulus;
     end 

 Vold=Vcalc;
 % calculate Iion using vector math now
 % You fill this part in
Ik= (gk*(ngate(1)^4)*(Vcalc-Ek));
Ina= (gna*(mgate(1)^3)*hgate(1)*(Vcalc-Ena));
Il= gl*(Vcalc-El);
Iion= Ik + Ina + Il;
 %Vector math and transcribing
 
 Vcalc=Vold(:) + (dt/Cm).*(coeff*A*Vold(:)-Iion(:)+ Istim(:));
 Vcalc=Vcalc.';
 smallVm=Vcalc-Vrest;

 outvm1(i)=Vcalc(1);
 outvm2(i)=Vcalc(50);
 outvm3(i)=Vcalc(100);

  mgate=mgate+dt.*(alpham(smallVm).*(1-mgate)-betam(smallVm).*mgate);
  hgate=hgate+dt.*(alphah(smallVm).*(1-hgate)-betah(smallVm).*hgate);
  ngate=ngate+dt.*(alphan(smallVm).*(1-ngate)-betan(smallVm).*ngate);
  
  %this section is to determine when to dump the output
   cnt=cnt+1; 
   figure(1000)  
   % Plot every 20 steps 
   if cnt==20
     hold off
     hh=subplot(1,1,1);
     plot(x,Vcalc), axis([0 x(end) -80 60]);
     set(get(hh,'XLabel'),'String','Space (cm)');
     set(get(hh,'YLabel'),'String','Voltage (mV)');
     timer=num2str(i*dt);
     set(get(hh,'Title'),'String',[timer,'  ','ms'] );
     cnt=0;
     drawnow % you need this to get MATLAB to draw each time
   end
   if(tvect(i) > time1 && pcount==0)
       Vprint1=Vcalc;
       pcount = 1;
   elseif(tvect(i)>time2 && pcount == 1) 
       Vprint2 = Vcalc;
       pcount = 2;
   elseif(tvect(i)>time3 && pcount == 2)
       Vprint3=Vcalc;
       pcount=3;
   end

 end

 figure;
 plot(x,Vprint1,'LineWidth',3), axis([0 x(end) -80 60]);
 hold on
 plot(x,Vprint2,'LineWidth',3)
 hold on
 plot(x,Vprint3,'LineWidth',3)
 
figure;
plot(tvect,outvm1,'LineWidth',3)
xlabel('time (msec)')
ylabel('Vm (mV)')
hold on
plot(tvect,outvm2,'LineWidth',3)
hold on
plot(tvect,outvm3,'LineWidth',3)
hold on

%end
