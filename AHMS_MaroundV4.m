clear all;


MyMat = load('ae4316p_2022_data_group3.mat')
% You can uncomment the functions one at the time. The matrix that needs to
% be entered in SPSS will be a 6x3 matrix and will appear in the Command
% window if you scroll up a bit after the function has run.

% Errorbar plot for dependent variable #1, Score Parameter
% DisplaySc(MyMat,2)

% Errorbar plot for dependent variable #2, Control signal variance
% DisplaySigU(MyMat,1)

%_________________ from this part, please modify the objesctivefcn program
%the for loop to run from 4 to iterator, not iterator -1.

% Errorbar plot for dependent variable #3, Pilot gain
% DisplayKp(MyMat,2)


% Errorbar plot for dependent variable #4, Pilot lead time constant - NOT
% NORMALLY distributed after performing the tests
% DisplayTL(MyMat,2)

% Errorbar plot for dependent variable #5, Pilot effective time delay
% DisplaytP (MyMat,2)

% Errorbar plot for dependent variable #6, Neuro-muscular natural frequency
% Displaywnm (MyMat,2)

% Errorbar plot for dependent variable #7, Neuro-muscular damping ratio
% Displayzetanm (MyMat,2)

% Errorbar plot for dependent variable #8, Crossover frequency
% DisplayWc (MyMat,2)

% Errorbar plot for dependent variable #9, Phase margin
% DisplayPm (MyMat,2)




%---------------- Plotting If needed ------

% w = 0.1:0.01:12
% 
% 
% [a,b] = AverageData_new(MyMat,'CM','data_subj1')
% 
% mag = a
% phase = b
% 
% mimi = optimize_simple(MyMat,mag,phase)
% [Hp1,Hp2] = get_simple_pilot(w,mimi)
% 
% [Hp3,Hp4] = open_loop_simple(Hp1,Hp2,w)
% 
% [crossover, phasemargin] = CrosPh_simple(Hp3(w<6),Hp4(w<6),w)
% 
% b(end) = b(end) - 360

hold on
plot(w,Hp1)
scatter(MyMat.data_subj1.C.w,a)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\omega [{rad/s}]')
ylabel('|H_{p}(j\omega)| [-]')
grid on
hold off

% 
% hold on
% scatter(MyMat.data_subj1.C.w,b)
% plot(w,Hp2)
% set(gca,'xscale','log')
% xlabel('\omega [{rad/s}]')
% ylabel('|H_{p}(j\omega)| [-]')
% hold off






function [crossover_freq,phase_margin] = CrosPh_simple(mags,angles,w) 


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


end





function elements = optimize_simple(data,mag,phase)


fun = @(params)objectivefcn(params,data.data_subj1.C.w,data,mag,phase);

params0 = [3,1.5,0.3,10,0.5];

options = optimset('MaxFunEvals',10000,'MaxIter',10000);
elements = fminsearch(fun,params0,options);


end

function [out1, out2] = get_simple_pilot(w,elements)

out = elements(1).*(elements(2).*w.*1j+1).*exp(-elements(3).*1j.*w).*(elements(4).^2)./((1j.*w).^2+2.*elements(4).*elements(5)*1j.*w+elements(4).^2);

out1 = abs(out)

out2 = phase(out)

ind = find(out2(end:-1:1) <0)

out2((numel(out2)-ind(1))+2:end) = out2((numel(out2)-ind(1))+2:end) - 2*pi

out2 = out2*180/pi


end

function [mag,phase] = open_loop_simple(magP,phaseP,omega)
    [magC, phaseC] = ControlledDyn(omega)
    mag = magP.*magC
    phase = phaseP + phaseC - 360
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
    hold on
    means = mean(MatSigSc)
    plot(MatSigSc','b-.x')
    
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('S_{c} [-]')
    hold off
else


    means = mean(MatSigSc)
    stdev = std(MatSigSc)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('S_{c} [-]')
    xlim([-0.3 4])
    ylim([0 5])
    grid on;
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
    hold on
    means = mean(MatSigU)
    plot(MatSigU','k-.s','MarkerFaceColor','black')
    plot(means,'r-s','MarkerFaceColor','red')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\sigma_{u}^2 [deg^2]')
    grid on
    hold off

else

    means = mean(MatSigU)
    stdev = std(MatSigU)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\sigma^2_{u} [deg^2]')
    xlim([-0.3 4])
    ylim([-10 80])
    grid on

end
end



function [] = DisplayKp(Data, choice)
names = fieldnames(Data);

%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6;
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = (names{i});

    [aC,bC] = AverageData_new(Data,'C',subj)
    mimiC = optimize_simple(Data,aC,bC)

    [aCP,bCP] = AverageData_new(Data,'CP',subj)
    mimiCP = optimize_simple(Data,aCP,bCP)

    [aCM,bCM] = AverageData_new(Data,'CM',subj)
    mimiCM = optimize_simple(Data,aCM,bCM)

    Cs = [Cs ; mimiC(1)];
    CPs = [CPs ; mimiCP(1)];
    CMs = [CMs ; mimiCM(1)];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

MatKp = [Cs,CPs,CMs]

display(MatKp)

if (choice == 1)
    plot(MatKp','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\sigma_{u}^2 [deg^2]')

else

    means = mean(MatKp)
    stdev =  std(MatKp)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('K_{p} [deg/deg]')
    xlim([-0.3 4])
    ylim([0 2])
    grid on;

end
end


function [] = DisplayTL(Data, choice)
names = fieldnames(Data);

%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6;
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = (names{i});

    [aC,bC] = AverageData_new(Data,'C',subj)
    mimiC = optimize_simple(Data,aC,bC)

    [aCP,bCP] = AverageData_new(Data,'CP',subj)
    mimiCP = optimize_simple(Data,aCP,bCP)

    [aCM,bCM] = AverageData_new(Data,'CM',subj)
    mimiCM = optimize_simple(Data,aCM,bCM)

    Cs = [Cs ; mimiC(2)];
    CPs = [CPs ; mimiCP(2)];
    CMs = [CMs ; mimiCM(2)];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

MatTL = [Cs,CPs,CMs]

display(MatTL)

if (choice == 1)
    plot(MatTL','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('T_{L} [sec]')

else

    means = mean(MatTL)
    stdev =  std(MatTL)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('T_{L} [sec]')
    xlim([-0.3 4])

end
end


function [] = DisplaytP(Data, choice)
names = fieldnames(Data);

%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6;
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = (names{i});

    [aC,bC] = AverageData_new(Data,'C',subj)
    mimiC = optimize_simple(Data,aC,bC)

    [aCP,bCP] = AverageData_new(Data,'CP',subj)
    mimiCP = optimize_simple(Data,aCP,bCP)

    [aCM,bCM] = AverageData_new(Data,'CM',subj)
    mimiCM = optimize_simple(Data,aCM,bCM)

    Cs = [Cs ; mimiC(3)];
    CPs = [CPs ; mimiCP(3)];
    CMs = [CMs ; mimiCM(3)];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

MattP = [Cs,CPs,CMs]

display(MattP)

if (choice == 1)
    plot(MattP','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\sigma_{u}^2 [deg^2]')

else

    means = mean(MattP)
    stdev =  std(MattP)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\tau_p [sec]')
    xlim([-0.3 4])
    ylim([0 1])
    grid on;

end
end

function [] = Displaywnm(Data, choice)
names = fieldnames(Data);

%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6;
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = (names{i});

    [aC,bC] = AverageData_new(Data,'C',subj)
    mimiC = optimize_simple(Data,aC,bC)

    [aCP,bCP] = AverageData_new(Data,'CP',subj)
    mimiCP = optimize_simple(Data,aCP,bCP)

    [aCM,bCM] = AverageData_new(Data,'CM',subj)
    mimiCM = optimize_simple(Data,aCM,bCM)

    Cs = [Cs ; mimiC(4)];
    CPs = [CPs ; mimiCP(4)];
    CMs = [CMs ; mimiCM(4)];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

Matwnm = [Cs,CPs,CMs]

display(Matwnm)

if (choice == 1)
    plot(MatSigU','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\sigma_{u}^2 [deg^2]')

else

    means = mean(Matwnm)
    stdev =  std(Matwnm)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\omega_{nm} [rad/s]')
    xlim([-0.3 4])
    ylim([5 12])
    grid on;

end
end

function [] = Displayzetanm(Data, choice)
names = fieldnames(Data);

%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6;
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = (names{i});

    [aC,bC] = AverageData_new(Data,'C',subj)
    mimiC = optimize_simple(Data,aC,bC)

    [aCP,bCP] = AverageData_new(Data,'CP',subj)
    mimiCP = optimize_simple(Data,aCP,bCP)

    [aCM,bCM] = AverageData_new(Data,'CM',subj)
    mimiCM = optimize_simple(Data,aCM,bCM)

    Cs = [Cs ; mimiC(5)];
    CPs = [CPs ; mimiCP(5)];
    CMs = [CMs ; mimiCM(5)];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

Matzetanm = [Cs,CPs,CMs]

display(Matzetanm)

if (choice == 1)
    plot(Matzetanm','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\sigma_{u}^2 [deg^2]')

else

    means = mean(Matzetanm)
    stdev =  std(Matzetanm)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\zeta_{nm} [-]')
    xlim([-0.3 4])
    ylim([0 1])
    grid on;

end
end



function [] = DisplayWc(Data, choice)
names = fieldnames(Data);
w = 0.1:0.001:12;

%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = Data.(names{i});

    arrC = []
    arrCP = []
    arrCM = []

    for j = 1:5

        optimizesC = optimize_simple(Data,subj.C.mag_Hp(:,j),subj.C.phase_Hp(:,j))
        optimizesCP = optimize_simple(Data,subj.CP.mag_Hp(:,j),subj.CP.phase_Hp(:,j))
        optimizesCM = optimize_simple(Data,subj.CM.mag_Hp(:,j),subj.CM.phase_Hp(:,j))

        [Hp1C,Hp2C] = get_simple_pilot(w,optimizesC)
        [Hp3C,Hp4C] = open_loop_simple(Hp1C,Hp2C,w)
        [crossoverC, phasemarginC] = CrosPh_simple(Hp3C(w<6),Hp4C(w<6),w)

        [Hp1CP,Hp2CP] = get_simple_pilot(w,optimizesCP)
        [Hp3CP,Hp4CP] = open_loop_simple(Hp1CP,Hp2CP,w)
        [crossoverCP, phasemarginCP] = CrosPh_simple(Hp3CP(w<6),Hp4CP(w<6),w)
        
        [Hp1CM,Hp2CM] = get_simple_pilot(w,optimizesCM)
        [Hp3CM,Hp4CM] = open_loop_simple(Hp1CM,Hp2CM,w)
        [crossoverCM, phasemarginCM] = CrosPh_simple(Hp3CM(w<6),Hp4CM(w<6),w)

        arrC = [arrC,crossoverC]
        arrCP = [arrC,crossoverCP]
        arrCM = [arrC,crossoverCM]

    end

    Cs = [Cs ; mean(arrC)];
    CPs = [CPs ; mean(arrCP)];
    CMs = [CMs ; mean(arrCM)];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

Matwc = [Cs,CPs,CMs]

display(Matwc)

if (choice == 1)
    plot(Matwc','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\omega_{c} [rad/s]')
    grid on

else

    means = mean(Matwc)
    stdev = std(Matwc)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\omega_{c} [rad/s]')
    xlim([-0.3 4])
    ylim([0 6])
    grid on

end
end


function [] = DisplayPm(Data, choice)
names = fieldnames(Data);
w = 0.1:0.001:12;

%loop through the subjects
Cs = [];
CPs = [];
CMs = [];

for i = 1:6
    % Append the mean of the variances for each subject and for each
    % condition C, CP, CM
    subj = Data.(names{i});

    arrC = []
    arrCP = []
    arrCM = []

    for j = 1:5

        optimizesC = optimize_simple(Data,subj.C.mag_Hp(:,j),subj.C.phase_Hp(:,j))
        optimizesCP = optimize_simple(Data,subj.CP.mag_Hp(:,j),subj.CP.phase_Hp(:,j))
        optimizesCM = optimize_simple(Data,subj.CM.mag_Hp(:,j),subj.CM.phase_Hp(:,j))

        [Hp1C,Hp2C] = get_simple_pilot(w,optimizesC)
        [Hp3C,Hp4C] = open_loop_simple(Hp1C,Hp2C,w)
        [crossoverC, phasemarginC] = CrosPh_simple(Hp3C(w<6),Hp4C(w<6),w)

        [Hp1CP,Hp2CP] = get_simple_pilot(w,optimizesCP)
        [Hp3CP,Hp4CP] = open_loop_simple(Hp1CP,Hp2CP,w)
        [crossoverCP, phasemarginCP] = CrosPh_simple(Hp3CP(w<6),Hp4CP(w<6),w)
        
        [Hp1CM,Hp2CM] = get_simple_pilot(w,optimizesCM)
        [Hp3CM,Hp4CM] = open_loop_simple(Hp1CM,Hp2CM,w)
        [crossoverCM, phasemarginCM] = CrosPh_simple(Hp3CM(w<6),Hp4CM(w<6),w)

        arrC = [arrC,phasemarginC]
        arrCP = [arrC,phasemarginCP]
        arrCM = [arrC,phasemarginCM]

    end

    Cs = [Cs ; mean(arrC)];
    CPs = [CPs ; mean(arrCP)];
    CMs = [CMs ; mean(arrCM)];

end

% Construct a final matrix where each line represents a subject, and column
% the different condition C, CP, CM

Matpm = [Cs,CPs,CMs]

display(Matpm)

if (choice == 1)
    plot(Matpm','b-x')
    nms = 1:1:3
    set(gca,'xTick',nms)
    set(gca,'xTickLabels',Data.str_conds)
    xlabel('Condition')
    ylabel('\phi_{m} [deg]')
    grid on

else

    means = mean(Matpm)
    stdev = std(Matpm)
    nms = 1:1:3
    errorbar(nms,means,stdev,'-s','MarkerFaceColor','black')
    set(gca,'xTick',nms)
    set(gca,'xTickLabels')
    xlabel('Condition')
    ylabel('\phi_{m} [deg]')
    xlim([-0.3 4])
    ylim([0 100])
    grid on

end
end




function [newmag,newphase] = AverageData_new(Data,choice,subj)

Mag = Data.(subj).(choice).mag_Hp
Phase = Data.(subj).(choice).phase_Hp

totaldata = Mag.*exp(1j.*deg2rad(Phase))


sumdata = sum(totaldata,2)
finaldata = sumdata/5


newmag = abs(finaldata)
newphase = rad2deg(angle(finaldata))


end



function [mag, phass] = ControlledDyn(omega)
Hc = 4./((omega*1j).^2)

mag = abs(Hc)
phass = angle(Hc) * 180/pi 


end


% function [newmag,newphase] = AverageData_new_v2(Data,choice,subj)
% 
% Mag = Data.(subj).(choice).mag_Hp
% Phase = Data.(subj).(choice).phase_Hp
% 
% totaldata = Mag.*exp(1j.*deg2rad(Phase))
% 
% 
% sumdata = sum(Mag,2)
% finaldata = sumdata/5
% 
% newmag = finaldata
% 
% sumdata2 = sum(Phase,2)
% finaldata2 = sumdata2/5
% 
% newphase = finaldata2
% end

% function [newmag,newphase,finaldata] = AverageData(Data,choice)
% 
% [Mag,Phase] = ConcatFresp(Data,choice)
% totaldata = Mag.*exp(1j.*deg2rad(Phase))
% 
% 
% sumdata = sum(totaldata,2)
% finaldata = sumdata/30
% 
% 
% newmag = abs(finaldata)
% newphase = rad2deg(angle(finaldata))
% 
% 
% end



% function [Mag, Phase] = ConcatFresp(Data,choice)
% names = fieldnames(Data);
% 
% Mag = []
% Phase = []
% for i = 1:6
%     Mag =  [Mag,Data.(names{i}).(choice).mag_Hp];
%     Phase = [Phase,Data.(names{i}).(choice).phase_Hp];
% 
% end
% end
% 


% function elements = optimize(data, choice)
% 
% [mag,phase] = AverageData(data,choice)
% 
% 
% fun = @(params)objectivefcn(params,data.data_subj1.C.w,data,mag,phase);
% 
% params0 = [3,1.5,0.3,10,0.5];
% 
% options = optimset('MaxFunEvals',10000,'MaxIter',10000);
% elements = fminsearch(fun,params0,options);
% 
% 
% end


% function dispCombined(w, data, whatplot)
% [abs_C,phase_C_fit, abs_C_pilot, phase_C_pilot] = plot_freq_domain(w,data,'C','Nothing')
% [abs_CP,phase_CP_fit, abs_CP_pilot, phase_CP_pilot] = plot_freq_domain(w,data,'CP','Nothing')
% [abs_CM,phase_CM_fit, abs_CM_pilot, phase_CM_pilot] = plot_freq_domain(w,data,'CM','Nothing')
% 
% [mag_C, phase_C] = AverageData(data,'C')
% [mag_CM, phase_CM] = AverageData(data,'CM')
% [mag_CP, phase_CP] = AverageData(data,'CP')
% 
% phase_C(end) = phase_C(end) - 360
% phase_CM(end) = phase_CM(end) - 360
% phase_CP(end) = phase_CP(end) - 360
% 
% if strcmp(whatplot,'Pilot')
% 
%     subplot(2,1,1)
%     hold on
%     plot(w, abs_C_pilot,'Color','k')
%     scatter(data.data_subj1.C.w, mag_C,'ko')
%     plot(w, abs_CP_pilot,"Color",'b')
%     scatter(data.data_subj1.C.w, mag_CP,'b^')
%     plot(w, abs_CM_pilot,"Color",'r')
%     scatter(data.data_subj1.C.w, mag_CM,'rd')
% 
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
%     xlabel('\omega [{rad/s}]')
%     ylabel('|H_{p}(j\omega)| [-]')
%     legend('C-fit','C-raw','CP-fit','CP-raw','CM-fit','CM-raw')
%     grid on
%     hold off
%     
%     
%     
%     subplot(2,1,2)
%     hold on
%     plot(w, phase_C_pilot,'Color','k')
%     scatter(data.data_subj1.C.w, phase_C,'ko')
%     plot(w, phase_CP_pilot,"Color",'b')
%     scatter(data.data_subj1.C.w, phase_CP,'b^')
%     plot(w, phase_CM_pilot,"Color",'r')
%     scatter(data.data_subj1.C.w, phase_CM,'rd')
%     set(gca,'xscale','log')
%     xlabel('\omega [{rad/s}]')
%     ylabel('\angle H_{p}(j\omega) [\circ]')
%     legend('C-fit','C-raw','CP-fit','CP-raw','CM-fit','CM-raw')
%     grid on
%     hold off
%     % 
% 
% elseif strcmp(whatplot,'Combined')
% 
% 
%     subplot(2,1,1)
%     hold on
%     plot(w,abs_C,w,abs_CP,w,abs_CM)
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
%     xlabel('\omega [{rad/s}]')
%     ylabel('|H_{p}(j\omega) \cdot H_{c}(j\omega)| [-]')
%     legend('Open loop - C', 'Open loop - CP', 'Open loop - CM')
%     grid on 
%     hold off
%     
%     
%     
%     subplot(2,1,2)
%     hold on
%     plot(w,phase_C_fit,w,phase_CP_fit,w,phase_CM_fit)
%     set(gca,'xscale','log')
%     xlabel('\omega [{rad/s}]')
%     ylabel('\angle H_{p}(j\omega) + H_{c}(j\omega) [\circ]')
%     legend('Open loop - C', 'Open loop - CP', 'Open loop - CM')
%     grid on
%     hold off
% 
% end
% end
% 
% function  dispPilotParams (data,choice)
% 
% params = optimize(data,choice)
% 
% display(['Pilot gain, Kp is : ', num2str(params(1))])
% display(['Lead time constant, TL is : ', num2str(params(2))])
% display(['Pilot time delay, tau_p is : ', num2str(params(3))])
% display(['Neuromuscular natural frequency, omega_nm is : ', num2str(params(4))])
% display(['Neuromuscular damping ratio, epsilon_nm is : ', num2str(params(5))])
% 
% end
% 
% function [total_abs, total_phase, pilot_abs, pilot_phase] = plot_freq_domain (w,data,choice,whatplot)
% 
% [mag,phase] = AverageData(data,choice)
% 
% pilot_params = optimize(data,choice)
% 
% yy = pilot_params(1).*(pilot_params(2).*w.*1j+1).*exp(-pilot_params(3).*1j.*w).*(pilot_params(4).^2)./((1j.*w).^2+2.*pilot_params(4).*pilot_params(5)*1j.*w+pilot_params(4).^2);
% 
% uu = angle(yy)
% 
% ind = find(uu(end:-1:1) <0)
% 
% uu((numel(uu)-ind(1))+2:end) = uu((numel(uu)-ind(1))+2:end) - 2*pi
% phase(end) = phase(end) - 360
% 
% [ss,nn] = ControlledDyn(w)
% 
% if strcmp(whatplot,'Pilot')
% 
%     subplot(2,1,1)
%     hold on
%     plot(w,abs(yy))
%     scatter(data.data_subj1.C.w,mag,'s','MarkerFaceColor',[0 0.447 0.741])
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
%     xlabel('\omega [{rad/s}]')
%     ylabel('|H_{p}(j\omega)| [-]')
%     legend('Fitted model','Raw data')
%     grid on
%     hold off
%     
%     
%     
%     subplot(2,1,2)
%     hold on
%     plot(w,uu.*180/pi)
%     scatter(data.data_subj1.C.w,phase,'s','MarkerFaceColor',[0 0.447 0.741])
%     set(gca,'xscale','log')
%     xlabel('\omega [{rad/s}]')
%     ylabel('\angle H_{p}(j\omega) [\circ]')
%     legend('Fitted model','Raw data')
%     grid on
%     hold off
%     % 
% 
% 
% elseif strcmp(whatplot,'Combined')
% 
% 
%     subplot(2,1,1)
%     hold on
%     plot(w,abs(yy).*ss)
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
%     xlabel('\omega [{rad/s}]')
%     ylabel('|H_{p}(j\omega) \cdot H_{c}(j\omega)| [-]')
%     grid on 
%     hold off
%     
%     
%     
%     subplot(2,1,2)
%     hold on
%     plot(w,uu.*180/pi+nn-360)
%     set(gca,'xscale','log')
%     xlabel('\omega [{rad/s}]')
%     ylabel('\angle H_{p}(j\omega) + H_{c}(j\omega) [\circ]')
%     grid on
%     hold off
% 
% end
% 
% total_abs = abs(yy).*ss
% total_phase = uu.*180/pi+nn-360
% 
% pilot_abs = abs(yy)
% pilot_phase = uu.*180/pi
% 
% 
% end
% function dispCrosPh (w,data,choice) 
% 
% [mags,angles,n,nn] = plot_freq_domain(w,data,choice,'Nothing')
% 
% idx1 = mags>1
% idx2 = mags<1
% 
% mag1 = mags(idx1)
% mag2 = mags(idx2)
% 
% w1 = w(idx1)
% w2 = w(idx2)
% 
% ang1 = angles(idx1)
% ang2 = angles(idx2)
% 
% mags_inter = [mag1(end), mag2(1)]
% w_inter = [w1(end), w2(1)]
% ang_inter = [ang1(end), ang2(1)]
% 
% wwi = 1
% crossover_freq = interp1(mags_inter,w_inter,wwi)
% 
% 
% 
% bb = interp1(w_inter,ang_inter,crossover_freq)
% 
% phase_margin = 180 + bb
% 
% display(['The crossover frequency is : ', num2str(crossover_freq),' [rad/s]'])
% display(['The phase margin is : ', num2str(phase_margin),' [deg]'])
% 
% end

