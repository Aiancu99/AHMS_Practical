clear all;
MyMat = load('ae4316p_2022_data_group3.mat')


% Find the name of the field variables of the struct

% DisplaySigU(MyMat,1)
% DisplaySc(MyMat,1)
% 

% 

% 
% w = 1:0.01:12;
% 
% yy = dd(1).*(dd(2).*w.*1j+1).*exp(-dd(3).*1j.*w).*(dd(4).^2)./((1j.*w).^2+2.*dd(4).*dd(5)*1j.*w+dd(4).^2);
% hold on
% plot(w,angle(yy)*180/pi)
% scatter(x,phase)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% 
% hold off

% scatter(MyMat.data_subj1.C.w,MyMat.data_subj2.C.mag_Hp(:,1))
% set(gca,'xscale','log')
% set(gca,'yscale','log')



% [mag,phase] = AverageData(MyMat,'CM')
% 
% 
% x=MyMat.data_subj2.C.w
% data = MyMat
% absol = mag
% phase = phase
% 
% fun = @(params)objectivefcn(params,x,data,absol,phase)
% 
% params0 = [3,1.5,0.3,10,0.5]
% 
% options = optimset('MaxFunEvals',4000,'MaxIter',4000);
% dd = fminsearch(fun,params0,options);
% 
% w = 0.1:0.01:12;
% yy = dd(1).*(dd(2).*w.*1j+1).*exp(-dd(3).*1j.*w).*(dd(4).^2)./((1j.*w).^2+2.*dd(4).*dd(5)*1j.*w+dd(4).^2);
% 
% subplot(2,1,1)
% hold on
% plot(w,abs(yy))
% scatter(MyMat.data_subj1.C.w,mag,'s','MarkerFaceColor',[0 0.447 0.741])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('\omega [{rad/s}]')
% ylabel('|H_{p}(j\omega)| [-]')
% legend('Fitted model','Raw data')
% hold off
% 
% 
% 
% subplot(2,1,2)
% hold on
% plot(w,angle(yy).*180/pi)
% scatter(MyMat.data_subj1.C.w,phase,'s','MarkerFaceColor',[0 0.447 0.741])
% set(gca,'xscale','log')
% xlabel('\omega [{rad/s}]')
% ylabel('\angle H_{p}(j\omega) [\circ]')
% legend('Fitted model','Raw data')
% hold off
% 

DisplaySc(MyMat,0)

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
    plot(MatSigU','b-x')
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


function [newmag,newphase] = AverageData(Data,choice)

[Mag,Phase] = ConcatFresp(Data,choice)
totaldata = Mag.*exp(1j.*deg2rad(Phase))
sumdata = sum(totaldata,2)
finaldata = sumdata/30

newmag = abs(finaldata)
newphase = rad2deg(phase(finaldata))


end