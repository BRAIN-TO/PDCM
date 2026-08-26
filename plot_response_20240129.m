
% region
k = 1;

% input index
u = 1;

% stimulus length
l1 = 10;
l2 = l1 * 16;


% event start time
s1 = [10, 30, 50, 100, 120, 140, 190, 210, 230, 280, 300, 320]; % correspond to x2
s2 = s1.*16; % correspond to x1

% x for plot
x1 = [DCM.U.dt:DCM.U.dt:length(DCM.U.u(:,1))*DCM.U.dt];
x2 = [1:DCM.v]*DCM.Y.dt;

% input
%plot(x1,DCM.U.u(:,u));
%hold on
% response
%plot(x2,DCM.xY(1).u(:,k));
%hold on
% response
%plot(x2,DCM.Y.y(:,k));
%hold on
% predicted
%plot(x2,DCM.y(:,k));
%hold on
% neuronal activity
%plot(x1, DCM.X(:,1));
%hold on
% rCBF
%plot(x1, DCM.X(:,3));
%hold on

for i = 1:length(s1)
    t1 = s1(i);
    t2 = s2(i);
    new_y(:,i) = DCM.y(t1-l1/2:t1+2*l1-1,k);
    new_Y(:,i) = DCM.Y.y(t1-l1/2:t1+2*l1-1,k);
    new_X1(:,i) = DCM.X(t2-l2/2:t2+2*l2-1,(k-1)*6+1);
    new_X3(:,i) = DCM.X(t2-l2/2:t2+2*l2-1,(k-1)*6+3);
end

figure('Position', [100 0 800 600]);
new_X1 = mean(new_X1,2);
plot(x1(1:400),new_X1(:,1).*10);
hold on

new_X3 = mean(new_X3,2);
plot(x1(1:400),new_X3(:,1).*10);
hold on


new_y = mean(new_y, 2);
plot(x2(1:2*l1+5),new_y(:,1)+0.6);
hold on

new_Y = mean(new_Y, 2);
%new_Y = new_Y + 1;
plot(x2(1:2*l1+5),new_Y(:,1)+1);


%plot(x1(1:415),DCM.U.u(81:495,u));
%legend
legend("estimated neuronal", "estimated CBF", "predicted BOLD", "measured BOLD");


%{
x3 = [2.898:2.898:1159.2];
s3 = floor(s2./14.4); % correspond to x3
l3 = 10;

for i = 1:length(s1)
    t1 = s1(i);
    t3 = s3(i);
    new_y(:,i) = DCM.y(t1-l1/2:t1+2*l1-1,k);
    new_Y(:,i) = DCM.Y.y(t1-l1/2:t1+2*l1-1,k);
    new_X1(:,i) = X((k-1)*5+2, t3-l3/2:t3+2*l3-1);
    new_X3(:,i) = X((k-1)*5+4, t3-l3/2:t3+2*l3-1);
end

figure('Position', [100 0 800 600]);
new_X1 = mean(new_X1,2);
plot(x3(1:25),new_X1(:,1).*10);
hold on

new_X3 = mean(new_X3,2);
plot(x3(1:25),new_X3(:,1).*10);
hold on


new_y = mean(new_y, 2);
plot(x2(1:2*l1+5),new_y(:,1));
hold on

new_Y = mean(new_Y, 2);
%new_Y = new_Y + 1;
plot(x2(1:2*l1+5),new_Y(:,1)+1);


%plot(x1(1:415),DCM.U.u(81:495,u));
%legend
legend("estimated neuronal", "estimated CBF", "predicted BOLD", "measured BOLD");
%}