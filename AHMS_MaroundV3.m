clear all;


MyMat = load('ae4316p_2022_data_group3.mat')

w = 0.1:0.001:12;

data = MyMat

choice = 'C'

choice2 = 'Pilot'

% optimize(data,choice)
% plot_freq_domain(w,data,choice,choice2)

%dispCombined(w,data,choice2)

%dispPilotParams(data,choice)

dispCrosPh(w,data,choice)



function dispCrosPh (w,data,choice) 

[mags,angles,n,nn] = plot_freq_domain(w,data,choice,'Nothing')

idx1 = mags>1
idx2 = mags<1

mag1 = mags(idx1)
mag2 = mags(idx2)

w1 = w(idx1)
w2 = w(idx2)

ang1 = angles(idx1)
ang2 = angles(idx2)

mags_inter = [mag1(end), mag2(1)]
w_inter = [w1(end), w2(1)]
ang_inter = [ang1(end), ang2(1)]

wwi = 1
crossover_freq = interp1(mags_inter,w_inter,wwi)



bb = interp1(w_inter,ang_inter,crossover_freq)

phase_margin = 180 + bb

display(['The crossover frequency is : ', num2str(crossover_freq),' [rad/s]'])
display(['The phase margin is : ', num2str(phase_margin),' [deg]'])

end

function dispCombined(w, data, whatplot)
[abs_C,phase_C_fit, abs_C_pilot, phase_C_pilot] = plot_freq_domain(w,data,'C','Nothing')
[abs_CP,phase_CP_fit, abs_CP_pilot, phase_CP_pilot] = plot_freq_domain(w,data,'CP','Nothing')
[abs_CM,phase_CM_fit, abs_CM_pilot, phase_CM_pilot] = plot_freq_domain(w,data,'CM','Nothing')

[mag_C, phase_C] = AverageData(data,'C')
[mag_CM, phase_CM] = AverageData(data,'CM')
[mag_CP, phase_CP] = AverageData(data,'CP')

phase_C(end) = phase_C(end) - 360
phase_CM(end) = phase_CM(end) - 360
phase_CP(end) = phase_CP(end) - 360

if strcmp(whatplot,'Pilot')

    subplot(2,1,1)
    hold on
    plot(w, abs_C_pilot,'Color','k')
    scatter(data.data_subj1.C.w, mag_C,'ko')
    plot(w, abs_CP_pilot,"Color",'b')
    scatter(data.data_subj1.C.w, mag_CP,'b^')
    plot(w, abs_CM_pilot,"Color",'r')
    scatter(data.data_subj1.C.w, mag_CM,'rd')

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('\omega [{rad/s}]')
    ylabel('|H_{p}(j\omega)| [-]')
    legend('C-fit','C-raw','CP-fit','CP-raw','CM-fit','CM-raw')
    grid on
    hold off
    
    
    
    subplot(2,1,2)
    hold on
    plot(w, phase_C_pilot,'Color','k')
    scatter(data.data_subj1.C.w, phase_C,'ko')
    plot(w, phase_CP_pilot,"Color",'b')
    scatter(data.data_subj1.C.w, phase_CP,'b^')
    plot(w, phase_CM_pilot,"Color",'r')
    scatter(data.data_subj1.C.w, phase_CM,'rd')
    set(gca,'xscale','log')
    xlabel('\omega [{rad/s}]')
    ylabel('\angle H_{p}(j\omega) [\circ]')
    legend('C-fit','C-raw','CP-fit','CP-raw','CM-fit','CM-raw')
    grid on
    hold off
    % 

elseif strcmp(whatplot,'Combined')


    subplot(2,1,1)
    hold on
    plot(w,abs_C,w,abs_CP,w,abs_CM)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('\omega [{rad/s}]')
    ylabel('|H_{p}(j\omega) \cdot H_{c}(j\omega)| [-]')
    legend('Open loop - C', 'Open loop - CP', 'Open loop - CM')
    grid on 
    hold off
    
    
    
    subplot(2,1,2)
    hold on
    plot(w,phase_C_fit,w,phase_CP_fit,w,phase_CM_fit)
    set(gca,'xscale','log')
    xlabel('\omega [{rad/s}]')
    ylabel('\angle H_{p}(j\omega) + H_{c}(j\omega) [\circ]')
    legend('Open loop - C', 'Open loop - CP', 'Open loop - CM')
    grid on
    hold off

end
end

function  dispPilotParams (data,choice)

params = optimize(data,choice)

display(['Pilot gain, Kp is : ', num2str(params(1))])
display(['Lead time constant, TL is : ', num2str(params(2))])
display(['Pilot time delay, tau_p is : ', num2str(params(3))])
display(['Neuromuscular natural frequency, omega_nm is : ', num2str(params(4))])
display(['Neuromuscular damping ratio, epsilon_nm is : ', num2str(params(5))])

end

function [total_abs, total_phase, pilot_abs, pilot_phase] = plot_freq_domain (w,data,choice,whatplot)

[mag,phase] = AverageData(data,choice)

pilot_params = optimize(data,choice)

yy = pilot_params(1).*(pilot_params(2).*w.*1j+1).*exp(-pilot_params(3).*1j.*w).*(pilot_params(4).^2)./((1j.*w).^2+2.*pilot_params(4).*pilot_params(5)*1j.*w+pilot_params(4).^2);

uu = angle(yy)

ind = find(uu(end:-1:1) <0)

uu((numel(uu)-ind(1))+2:end) = uu((numel(uu)-ind(1))+2:end) - 2*pi
phase(end) = phase(end) - 360

[ss,nn] = ControlledDyn(w)
if strcmp(whatplot,'Pilot')

    subplot(2,1,1)
    hold on
    plot(w,abs(yy))
    scatter(data.data_subj1.C.w,mag,'s','MarkerFaceColor',[0 0.447 0.741])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('\omega [{rad/s}]')
    ylabel('|H_{p}(j\omega)| [-]')
    legend('Fitted model','Raw data')
    grid on
    hold off
    
    
    
    subplot(2,1,2)
    hold on
    plot(w,uu.*180/pi)
    scatter(data.data_subj1.C.w,phase,'s','MarkerFaceColor',[0 0.447 0.741])
    set(gca,'xscale','log')
    xlabel('\omega [{rad/s}]')
    ylabel('\angle H_{p}(j\omega) [\circ]')
    legend('Fitted model','Raw data')
    grid on
    hold off
    % 


elseif strcmp(whatplot,'Combined')


    subplot(2,1,1)
    hold on
    plot(w,abs(yy).*ss)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('\omega [{rad/s}]')
    ylabel('|H_{p}(j\omega) \cdot H_{c}(j\omega)| [-]')
    grid on 
    hold off
    
    
    
    subplot(2,1,2)
    hold on
    plot(w,uu.*180/pi+nn-360)
    set(gca,'xscale','log')
    xlabel('\omega [{rad/s}]')
    ylabel('\angle H_{p}(j\omega) + H_{c}(j\omega) [\circ]')
    grid on
    hold off

end

total_abs = abs(yy).*ss
total_phase = uu.*180/pi+nn-360

pilot_abs = abs(yy)
pilot_phase = uu.*180/pi


end




function elements = optimize(data, choice)

[mag,phase] = AverageData(data,choice)


fun = @(params)objectivefcn(params,data.data_subj1.C.w,data,mag,phase);

params0 = [3,1.5,0.3,10,0.5];

options = optimset('MaxFunEvals',10000,'MaxIter',10000);
elements = fminsearch(fun,params0,options);


end


function [] = DisplaySc(Data, choice)
names = fieldnames(Data);



%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6;
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = Data.(names{i});

    Cs = [Cs ; mean(var(subj.C.e))/var(subj.C.ft)];
    CPs = [CPs ; mean(var(subj.CP.e))/var(subj.CP.ft)];
    CMs = [CMs ; mean(var(subj.CM.e))/var(subj.CM.ft)];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

MatSigSc = [Cs,CPs,CMs]

if (choice == 1)
    plot(MatSigSc','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('S_{c} [-]')

else


    means = mean(MatSigSc)
    stdev = 2 .* std(MatSigSc)
    nms = 1:1:3
    errorbar(nms,means,stdev)
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('S_{c} [-]')
    xlim([-0.3 4])

end
end


function [] = DisplaySigU(Data, choice)
names = fieldnames(Data);

%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6;
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = Data.(names{i});

    Cs = [Cs ; mean(var(subj.C.u))];
    CPs = [CPs ; mean(var(subj.CP.u))];
    CMs = [CMs ; mean(var(subj.CM.u))];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

MatSigU = [Cs,CPs,CMs]

display(MatSigU)

if (choice == 1)
    plot(MatSigU','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\sigma_{u}^2 [deg^2]')

else

    means = mean(MatSigU)
    stdev = 2 .* std(MatSigU)
    nms = 1:1:3
    errorbar(nms,means,stdev)
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\sigma_{u} ^ 2 [deg^2]')
    xlim([-0.3 4])

end
end


function [Mag, Phase] = ConcatFresp(Data,choice)
names = fieldnames(Data);

Mag = []
Phase = []
for i = 1:6
    Mag =  [Mag,Data.(names{i}).(choice).mag_Hp];
    Phase = [Phase,Data.(names{i}).(choice).phase_Hp];

end
end


function [newmag,newphase,finaldata] = AverageData(Data,choice)

[Mag,Phase] = ConcatFresp(Data,choice)
totaldata = Mag.*exp(1j.*deg2rad(Phase))


sumdata = sum(totaldata,2)
finaldata = sumdata/30


newmag = abs(finaldata)
newphase = rad2deg(angle(finaldata))


end

function [mag, phass] = ControlledDyn(omega)
Hc = 4./((omega*1j).^2)

mag = abs(Hc)
phass = angle(Hc) * 180/pi

end

