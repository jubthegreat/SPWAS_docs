function FWdoubler(Rload, on_time, off_time, C1, C2)
%%%%%%%% ORIGINALLY WRITTEN BY PROFESSOR ANDREW HOLMES %%%%%%%%%%%%%
%Modified by Shaun Price to use individual PMOS and NMOS Rdson, and altered
%diode model to fit the RB531SM-30FH. The script was also made into a
%function which is to be called from a main script providing simulation
%parameters.

% models active full-wave voltage doubler with Schottky diodes + mosfets
% mosfet switching using monostables as in Harrasi circuit
% assumes sinusoidal source with source resistance Rs

% params
time_step=1e-6;     % nominal time step (S)
fs=500;     % source frequency (Hz)     NB 494 Hz gives 500 mVrms
Kgen=2.7e-4;    % rms source voltage per RPM (Intel prototype MWT) (V)
Vpk=sqrt(2)*Kgen*3.75*fs;      % source peak amplitude (V)
Rs=21;      % source impedance (Intel prototype MWT) (Ohms)
V0=0.026;   % effective thermal voltage of diode  (V)
I0=250e-6;    % reverse saturation current of diode  (A)
Rdsonn = 0.6; % N mosfet on-state resistance (Ohms)
Rdsonp = 1.5; % P mosfet on-state resistance (Ohms)

Voffset=0.03;  % mosfet switching comparator offset (V)

% calculate I-V pairs for (diode + Rs)
Imax=Vpk/Rs;    % max current in look-up table
Istep=1e-6;     % granularity in look-up table (current)
Idiode=0:Istep:Imax;
Vdiode=V0*log(1+Idiode/I0);    
IVtable=[Idiode ; Vdiode + Idiode*Rs];    % Schottky diodes   
%IVtable=[Idiode ; Idiode*Rs];   % ideal active rectifier
figure(1);
clf;
hold on;
plot(IVtable(2,:), IVtable(1,:));%Plot voltage vs current 
xlabel('Bias voltage');
ylabel('Diode current');

% set up loop parameters
total_cyles=100;
plot_cycles=5; 
samples_per_cycle=floor(1/(time_step*fs));
nsamp=samples_per_cycle*total_cyles;
nplot=plot_cycles*samples_per_cycle;
dt=1/(samples_per_cycle*fs)    % actual time step (factor of period)
t=0:nsamp-1;
t=dt*t;
Vs=Vpk*sin(2*pi*fs*t);
on_steps=floor(on_time/dt);     % mosfet on-time in steps
off_steps=floor(off_time/dt);   % mosfet off-time in steps

% initialisation
Is=zeros(1,nsamp); %source current (?)
Id1=zeros(1,nsamp); %diode1 current(?)
Id2=zeros(1,nsamp); %diode2 current(?)
Iout=zeros(1,nsamp); %Output current
Vc1=zeros(1,nsamp); %Voltage across capacitor 1(?)
Vc2=zeros(1,nsamp); %voltage across capacitor 2(?)
Vout=zeros(1,nsamp); %Output voltage
Q1=0; Q2=0;   % Both mosfets off initially
Q1_timer=off_steps; Q2_timer=off_steps;

for i=1:nsamp-1
% NB: calculation assumes non-overlapping diode currents
Iout(i)=Vout(i)/Rload;          % RESISTIVE LOAD
%Iout(i)=2.6e-3;                 % CONST CURRENT LOAD
if Q1==0    
    if Vs(i)>Vc1(i) %When source voltage greater than capacitor voltage but mosfet is off, diode conducts slightly
        [Y I]=min(abs(IVtable(2,:)-(Vs(i)-Vc1(i)))); %returns value and index of closest value in table
        Id1(i)=IVtable(1,I); %branch current into capacitor through diode
    end
    if Q1_timer < off_steps
        Q1_timer=Q1_timer+1;        % mosfet dead time
    else
        if Vs(i)>Vc1(i)+Voffset
            Q1=1;
            Q1_timer=0;
        end
    end   
else
    Id1(i)=(Vs(i)-Vc1(i))/(Rs+Rdsonn); %mosfet is on, so simple ohms law over Rds
    if Q1_timer < on_steps
        Q1_timer=Q1_timer+1;
    else
        Q1=0;
        Q1_timer=0;
    end
end
if Q2==0
    if Vs(i)<-Vc2(i)
        [Y I]=min(abs(IVtable(2,:)-(-Vc2(i)-Vs(i))));
        Id2(i)=IVtable(1,I);
    end
    if Q2_timer < off_steps
        Q2_timer=Q2_timer+1;        % mosfet dead time
    else
        if Vs(i)<-Vc2(i)-Voffset
            Q2=1;
            Q2_timer=0;
        end
    end
else
    Id2(i)=(-Vc2(i)-Vs(i))/(Rs+Rdsonp);
    if Q2_timer < on_steps
        Q2_timer=Q2_timer+1;
    else
        Q2=0;
        Q2_timer=0;
    end
end  
Is(i)=Id1(i)-Id2(i); %source current is sum of each branch current
Vc1(i+1)=Vc1(i)+(Id1(i)-Iout(i))*dt/C1; %Assume constant current into capacitor over short time interval
Vc2(i+1)=Vc2(i)+(Id2(i)-Iout(i))*dt/C2; %Assume constant current into capacitor over short time interval
Vout(i+1)=Vc1(i+1)+Vc2(i+1); %Vout is sum of capacitor voltages
end
Vgen=Vs-Is*Rs;

t_plot=t(nsamp-nplot+1:nsamp);
Vs_plot=Vs(nsamp-nplot+1:nsamp);
Vc1_plot=Vc1(nsamp-nplot+1:nsamp);
Vc2_plot=Vc2(nsamp-nplot+1:nsamp);
Is_plot=Is(nsamp-nplot+1:nsamp);
Pin_plot=Vs_plot.*Is_plot;
Vout_plot=Vout(nsamp-nplot+1:nsamp);
Iout_plot=Iout(nsamp-nplot+1:nsamp);
Pout_plot=Vout_plot.*Iout_plot;
Vgen_plot=Vgen(nsamp-nplot+1:nsamp);
Iout_average = mean(Iout)

figure(2);
clf;
subplot(3,1,1);
hold on;
plot(t_plot,Vout_plot,'r');
plot(t_plot,Vc1_plot,'g');
plot(t_plot,-Vc2_plot,'c');
legend('Vout', 'V Upper Capacitor', 'V Lower Capacitor');
%plot(t_plot,Vgen_plot,'m--');

subplot(3,1,2);
hold on;
plot(t_plot,Vs_plot,'b'); %Source voltage
plot(t_plot,Vgen_plot,'m'); %Voltage across 

legend('Vturbine', 'terminal voltage'); 
subplot(3,1,3);
hold on;
plot(t_plot,Pin_plot,'b');
plot(t_plot,Pout_plot,'r');
legend('Input power', 'Output power')

Pin_avg=sum(Pin_plot)/size(Pin_plot,2)  
Pout_avg=sum(Pout_plot)/size(Pout_plot,2)
Efficiency=100*Pout_avg/Pin_avg
Vout_min=min(Vout_plot)
end
    
    



