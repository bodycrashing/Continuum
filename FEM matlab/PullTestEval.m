%% Poissons forhold calculator
close all; clear; clc
for m=1:5
    p=[0 5 6 7 8];
    p=num2str(p(m));
    %Choose porosity
    hulst=(p);
    %Unit cell size
    uc=10;
    %Loops over pull tests (3 for each kind of specimen)
    for j=1:3
        %% Finds how many samples were conducted at given pull-test
        testnr = num2str(j);
        testfolder = strcat('C:\Maskin\9. semester\Continuum Mechanics\Experiment\data\DIC\',hulst,'mm',testnr);
        dat = dir(testfolder);
        testamount = size(dat,1)-3;
        
        %% Finds appropriate sampling points
        i = testamount;
        Tstamp = num2str(i,'%05i');
        file = strcat(hulst,'mm',testnr,'-',Tstamp,'_1');
        importfile(file);
        if str2num(hulst)~= 0
            [hz,vl] = samplingpoint(eyy,hulst);
        else
            hz=[30 70];
            vl=[40 160];
        end
        
        xdist(j,m) = abs(x_c(1,hz(2))-x_c(1,hz(1)));
        ydist(j,m) = abs(y_c(vl(2),1)-y_c(vl(1),1));
        
        % Loops over the time samples of the specific pull test
        for i=0:testamount
            % Imports data
            Tstamp=num2str(i,'%05i');
            file=strcat(hulst,'mm',testnr,'-',Tstamp,'_1');
            importfile(file);
            
            % Data treatment
            %Strain in the pull direction
            eyy_c = abs(v_c(vl(1), hz(1)) - v_c(vl(2),hz(1)) + v_c(vl(1),hz(2)) - v_c(vl(2),hz(2)))/(2*ydist(j,m));
            %Strain in the transverse direction
            exx_c = abs(u_c(vl(1), hz(2)) - u_c(vl(1),hz(1)) + u_c(vl(2), hz(2)) - u_c(vl(2),hz(1)))/(2*xdist(j,m));
            
            %Poisson Ration
            nu = exx_c/eyy_c;
            
            mnu(i+1,j+m*3-3)=nu;
            if i==testamount
                E(j,m)=force(j,m)/(30*3)/eyy_c;
            end
            
            %specificplots(i,j,v_c)
        end
        
        [nuaver(j,m)]=poiss(mnu(:,j+m*3-3));
    end
end

nu=mean(nuaver);
standrNU=std(nuaver);
dens=[0/10^2,2.5^2*pi/10^2,3^2*pi/10^2,3.5^2*pi/10^2,4^2*pi/10^2];
mnu
nu
standrE=std(E)
E=mean(E)
resultplot(dens,nu,E,standrE,mnu,standrNU)

function [x,y] = samplingpoint(eyy,hulst)
%Finds the sampling points appropriate considering the unit cell
%Smoothing out in the y-direction
windowSize = 5;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
eyy_s=filter(b,a,eyy');
eyy_s=eyy_s';

%subplot(1,2,1)
splot(eyy_s)

%Setting threshold as 40% of the maximum array value, though with max
%value of 0.02 to combat outliers.
threshold=min([max(max(eyy_s))/2.5,0.015]);
test=eyy_s>threshold;
[centers,radii,metric]=imfindcircles(test,[5 15]);
lr=length(radii);

%Finding corner holes by taking a mean of circle center values close to
%the most extreme circle positions.
y_min=mean(centers(centers(:,2)<min(centers(:,2))+15,2));
y_max=mean(centers(centers(:,2)>max(centers(:,2))-15,2));
if y_max-y_min<45
    y_min=y_min-(y_max-y_min);
end
%nasty special code to fix Ø8 holes
if hulst=='8' && lr==5
    y_min=35;
else if hulst=='8' && lr==8
        y_max=157;
    end
end

y_off=1/8*(y_max-y_min);

x_min=mean(centers(centers(:,1)<min(centers(:,1))+15,1));
x_max=mean(centers(centers(:,1)>max(centers(:,1))-15,1));

if x_max-x_min<45
    x_min=x_min-(x_max-x_min);
end

%nasty special code to fix Ø8 holes
if hulst=='8' && lr==5
    x_min=15;
    x_max=84;
else if hulst=='8' && lr==8
        x_min=15;
        x_max=78;
    end
end


x_off=1/4*(x_max-x_min);

x=round([x_min+x_off,x_max-x_off]);
y=round([y_min+y_off,y_max-y_off]);

%Plotting the found circles
%subplot(1,2,2)
imshow(test)
%title('example circle finder')

if lr>15;
    lr=15;
end
centersStrong15 = centers(1:lr,:);
radiiStrong15 = radii(1:lr);
metricStrong15 = metric(1:lr);
viscircles(centersStrong15, radiiStrong15,'EdgeColor','b');

%subplot(1,2,2)
hold on
plot(x,y,'ro')
plot(fliplr(x),y,'ro')
hold off

end

function splot(ImpVar)
%Plots a surface plot of desired variable
x = 1:size(ImpVar,2);
y = 1:size(ImpVar,1);
[X,Y] = meshgrid(x,y);
surf(x,y,ImpVar)
colorbar
end

function specificplots(i,j,v_c)
subplot(3,3,i+1)
splot(v_c)
hold on
p=num2str(i);
k=num2str(j);
test=strcat('timeinstant-',p,'-testsample-',k);
title(test)
view(90,0)
end

function resultplot(dens,nuaver,E,standrE,mnu,standrNu)
%Plotting overview of E and Nu
f1=figure;
subplot(2,1,1);
hold on
plot(dens,E,'g*');
plot(dens,E+standrE,'r*')
plot(dens,E-standrE,'r*')

% pol=polyfit(dens,E,2);
% ypol=polyval(pol,dens);
% plot(dens,ypol);


[f,g]=fit(dens',E','poly1')
plot(f,dens',E')
osti='f(x) = ax+b'
text(0.1,0.18,osti,'Units','normalized')
ci=confint(f)
koef=['a';'b']
koefval=num2str([f.p1;f.p2])
osti2=strcat(koef,' = ',koefval,'        ',' conf int 95%',' = ',num2str(ci'))
text(0.25,0.18,osti2,'Units','normalized')
title('Modulus of elasticity as E-bar')
xlabel('Porosity')
ylabel('Modulus of elasticity [GPa]')
axis([0 0.6 0 3])

subplot(2,1,2)
hold on
plot(dens,nuaver,'g*')
plot(dens,nuaver+standrNu,'r*')
plot(dens,nuaver-standrNu,'r*')
[f,g]=fit(dens',nuaver','poly1')
plot(f,dens',nuaver')
osti='f(x) = ax+b'
text(0.1,0.18,osti,'Units','normalized')
ci=confint(f)
koef=['a';'b']
koefval=num2str([f.p1;f.p2])
osti2=strcat(koef,' = ',koefval,'       ',' conf int 95%',' = ',num2str(ci'))
text(0.25,0.18,osti2,'Units','normalized')
title('Poissons ratio as nu-bar')
xlabel('Porosity')
ylabel('Poisson ratio')
axis([0 0.6 0 0.5])

%Scatterplot of Nu
f2=figure;
s=4
x=1:s+2;
mnu(mnu==0)=nan;
hold on
p = plot(x,mnu(s:end,1:3),':ro',x,mnu(s:end,4:6),':b*',x,mnu(s:end,7:9),':gd',x,mnu(s:end,10:12),':ms',x,mnu(s:end,13:15),':c+','LineWidth',1.1);
xlabel('timestep')
ylabel('poisson ratio')
axis([1 7 0 0.6])
legend('No holes 1','No holes 2','No holes 3','Ø5 1','Ø5 2','Ø5 3','Ø6 1','Ø6 2','Ø6 3','Ø7 1','Ø7 2','Ø7 3','Ø8 1','Ø8 2','Ø8 3')
end

function [nu]=poiss(filnu)
%Removing NaN and 0 values, taking mean and removing any outliers from
%that.

filnu(isnan(filnu))=[];
filnu(find(filnu==0))=[];
filnu=filnu(length(filnu)-3:end);
nu=mean(filnu);

end

function f=force(testnr,hulst)
%Maximum values of force from pull test, fits with final
ftabel=[0.7705 0.6849 0.6261 0.6122 0.4528;
    0.6689 0.6264 0.6382 0.5866 0.5378;
    0.6546 0.6260 0.5721 0.6168 0.5054];
f=ftabel(testnr,hulst);
end

