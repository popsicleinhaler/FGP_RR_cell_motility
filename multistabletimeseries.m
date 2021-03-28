% UZR={'delta_P'=0.04}

% global delta_P
% delta_P = 0.0004;
% global B
% B = 2;

spc = ["PaxFAK","Paxs","PaxsFAK","PaxGIT","PaxsGIT","Rac","Rho"];

%from 0
one = 1;
two = 2;
if two>one
    figure;
    tiledlayout(1,3);
    nexttile;
    [t,y1] = ode15s(@rhs_fun,[0,1],[0,0,0,0,0,0,0]);
    plot(t,y1,'linewidth',5);
    legend('PaxFAK0','Paxs','PaxsFAK','PaxGIT','PaxsGIT','Rac','Rho');
    title('from 0 transient');
    nexttile;
    [t,y80] = ode15s(@rhs_fun,[0,80],[0,0,0,0,0,0,0]);
    plot(t,y80,'linewidth',5);
    legend('PaxFAK0','Paxs','PaxsFAK','PaxGIT','PaxsGIT','Rac','Rho');
    title('less transient');
    nexttile;
    [t,y400] = ode15s(@rhs_fun,[0,400],[0,0,0,0,0,0,0]);
    plot(t,y400,'linewidth',5);
    legend('PaxFAK0','Paxs','PaxsFAK','PaxGIT','PaxsGIT','Rac','Rho');
    title('stable');



    %from 1
    figure;
    tiledlayout(1,3);
    nexttile;
    [t,y1] = ode15s(@rhs_fun,[0,1],[1,1,1,1,1,1,1]);
    plot(t,y1,'linewidth',5);
    legend('PaxFAK1','Paxs','PaxsFAK','PaxGIT','PaxsGIT','Rac','Rho');
    title('start from 1');
    nexttile;
    [t,y80] = ode15s(@rhs_fun,[0,80],[1,1,1,1,1,1,1]);
    plot(t,y80,'linewidth',5);
    legend('PaxFAK1','Paxs','PaxsFAK','PaxGIT','PaxsGIT','Rac','Rho');
    title('start from 1');
    nexttile;
    [t,y400] = ode15s(@rhs_fun,[0,400],[1,1,1,1,1,1,1]);
    plot(t,y400,'linewidth',5);
    legend('PaxFAK1','Paxs','PaxsFAK','PaxGIT','PaxsGIT','Rac','Rho');
    title('start from 1');
end

% for b = 2, delta_P=0.0004
UZ1=[0.24608629126, 0.19681006037, 0.1950735664, 0.0082636687772, 0.32753212942, 0.23448393114, 0.19901606098];
UZ2=[0.27335402958, 0.18770477847, 0.18598039315, 0.009248518871, 0.31461822208, 0.13321252216, 0.3742675298]; %unstable
UZ3=[0.31794055221, 0.17264533548, 0.17095721709, 0.010892682976, 0.29285077906, 0.088753898638, 0.45163935709];
UZ4=[0.67836711314, 0.02985781125, 0.029435142727, 0.026373891646, 0.057219700795, 0.074603753557, 0.48986094587];%unstable
UZ5=[0.69319774255, 0.00091911378556, 0.00090635457267, 0.027703302867, 0.0018111004761, 0.073765299752, 0.49749465465];

% [0.24608629126, 0.19681006037, 0.1950735664, 0.0082636687772, 0.32753212942, 0.23448393114, 0.19901606098]
% [0.27335402958, 0.18770477847, 0.18598039315, 0.009248518871, 0.31461822208, 0.13321252216, 0.3742675298]
% [0.31794055221, 0.17264533548, 0.17095721709, 0.010892682976, 0.29285077906, 0.088753898638, 0.45163935709]
% [0.67836711314, 0.02985781125, 0.029435142727, 0.026373891646, 0.057219700795, 0.074603753557, 0.48986094587]
% [0.69319774255, 0.00091911378556, 0.00090635457267, 0.027703302867, 0.0018111004761, 0.073765299752, 0.49749465465]

% for b =5, delta_P=0.06
% UZ5 = [0.69319773835, 0.000919108459, 0.00090634932023, 0.027703302842, 0.0018110899899, 0.07376529993, 0.4974946558];
% UZ4 = [0.67836539167, 0.029859043305, 0.029436357643, 0.026373794197, 0.057221996282, 0.074603834746, 0.48986059371];
% UZ3 = [0.2733514885, 0.1877056543, 0.18598126719, 0.009248426197, 0.31461947273, 0.13321223666, 0.37426795061];
% UZ2 = [0.074603751701, 0.48986095393, 0.029857783065, 0.029435114934, 0.057219648282, 0.67836715252, 0.026373893875];
% UZ1 = [0.24608631127, 0.19681005346, 0.19507355951, 0.008263669496, 0.3275321197, 0.23448393099, 0.19901606182];

span = 2500;

%start at stable states
figure;
hold on;
[t,y] = ode15s(@rhs_fun,[0,span],UZ1);
plot(t,y(:,6),'linewidth',5);
[t,y] = ode15s(@rhs_fun,[0,span],UZ3);
plot(t,y(:,6),'linewidth',5);
[t,y] = ode15s(@rhs_fun,[0,span],UZ5);
plot(t,y(:,6),'linewidth',5);

legend('top','middle','bottom')
hold off

%start at unstable states
figure;
hold on;
[t,y] = ode15s(@rhs_fun,[0,span],UZ2);
plot(t,y(:,6),'linewidth',5);
[t,y] = ode15s(@rhs_fun,[0,span],UZ4);
plot(t,y(:,6),'linewidth',5);

legend('top','bottom')
hold off

%start at some fraction of a vector toward stable state away from a saddle
f = 0.01;

twoone = UZ1 - UZ2;
twothree = UZ3 - UZ2;

IC1 = UZ2 + f*twoone;
IC3p = UZ2 + f*twothree;

fourthree = UZ3-UZ4;
fourfive = UZ5-UZ4;

IC3d = UZ4 + f*fourthree;
IC5 = UZ4 + f*fourfive;

[t1,y1] = ode15s(@rhs_fun,[0,span],IC1);
[t3p,y3p] = ode15s(@rhs_fun,[0,span],IC3p);
[t3d,y3d] = ode15s(@rhs_fun,[0,span],IC3d);
[t5,y5] = ode15s(@rhs_fun,[0,span],IC5);

figure
tiledlayout(3,3)
for i = 1:7
    nexttile
    hold on
    name = spc(i);
    plot(t1,y1(:,i),'linewidth',5);
    plot(t3p,y3p(:,i),'linewidth',5);
    plot(t3d,y3d(:,i),'linewidth',5);
    plot(t5,y5(:,i),'linewidth',5);
    legend('top','middle high','middle low','bottom')
    title(name)
    hold off
end
% epsilon = 0.1
% 
% i=0
% while i < 100
%     i=i+1
%     many = rnadi(1,7)
%     for i 1:many
%         which = randi(1,7)
%         
%         
%     
%     

