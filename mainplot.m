clear all;
close all;

%前5000个数据代表50ms
result11 = xlsread('ECN=5KB200KB10%,PFC=300K838/2F+20G+20G.xlsx');
result12 = xlsread('ECN=5KB200KB10%,PFC=300K838/2F+40G+40G.xlsx');
result13 = xlsread('ECN=5KB200KB10%,PFC=300K838/4F+10G+10G.xlsx');
result14 = xlsread('ECN=5KB200KB10%,PFC=300K838/4F+40G+40G.xlsx');
result15 = xlsread('ECN=5KB200KB10%,PFC=300K838/8F+5G+5G.xlsx');
result16 = xlsread('ECN=5KB200KB10%,PFC=300K838/8F+40G+40G.xlsx');

result21 = xlsread('ECN=5KB200KB10%,PFC=800K838/2F+20G+20G.xlsx');
result22 = xlsread('ECN=5KB200KB10%,PFC=800K838/2F+40G+40G.xlsx');
result23 = xlsread('ECN=5KB200KB10%,PFC=800K838/4F+10G+10G.xlsx');
result24 = xlsread('ECN=5KB200KB10%,PFC=800K838/4F+40G+40G.xlsx');
result25 = xlsread('ECN=5KB200KB10%,PFC=800K838/8F+5G+5G.xlsx');
result26 = xlsread('ECN=5KB200KB10%,PFC=800K838/8F+40G+40G.xlsx');

result31 = xlsread('ECN=200KB400KB10%,PFC=800K838/2F+20G+20G.xlsx');
result32 = xlsread('ECN=200KB400KB10%,PFC=800K838/2F+40G+40G.xlsx');
result33 = xlsread('ECN=200KB400KB10%,PFC=800K838/4F+10G+10G.xlsx');
result34 = xlsread('ECN=200KB400KB10%,PFC=800K838/4F+40G+40G.xlsx');
result35 = xlsread('ECN=200KB400KB10%,PFC=800K838/8F+5G+5G.xlsx');
result36 = xlsread('ECN=200KB400KB10%,PFC=800K838/8F+40G+40G.xlsx');

result41 = xlsread('ECN=200KB400KB50%,PFC=800K838/2F+20G+20G.xlsx');
result42 = xlsread('ECN=200KB400KB50%,PFC=800K838/2F+40G+40G.xlsx');
result43 = xlsread('ECN=200KB400KB50%,PFC=800K838/4F+10G+10G.xlsx');
result44 = xlsread('ECN=200KB400KB50%,PFC=800K838/4F+40G+40G.xlsx');
result45 = xlsread('ECN=200KB400KB50%,PFC=800K838/8F+5G+5G.xlsx');
result46 = xlsread('ECN=200KB400KB50%,PFC=800K838/8F+40G+40G.xlsx');


%画6行图，分别对应不同流数量、不同初始速度
%画4列图，每列对应吞吐，S0队列、S1队列、累积丢包

figure('color','w')
set(gcf,'Position',[100 100 3000 3000]);
set(gca,'position', [0.1 0.15 0.9 0.85]);

%单流吞吐
subplot(6,4,1)
plot(10:10:50000,result11(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result21(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result31(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result41(1:5000,1),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('个流吞吐 (Gbps)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 25])
set(gca,'FontName','Times New Roman','FontSize',14);

%S0队列
subplot(6,4,2)
plot(10:10:50000,result11(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result21(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result31(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result41(1:5000,2),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S0队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 500])
set(gca,'FontName','Times New Roman','FontSize',14);

%S1队列
subplot(6,4,3)
plot(10:10:50000,result11(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result21(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result31(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result41(1:5000,3),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S1队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 500])
set(gca,'FontName','Times New Roman','FontSize',14);

%丢包
subplot(6,4,4)
plot(10:10:50000,result11(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result21(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result31(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result41(1:5000,4),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('累积丢包 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 max(5,max(result11(1:5000,4)))])
set(gca,'FontName','Times New Roman','FontSize',14);



%第二行
%单流吞吐
subplot(6,4,5)
plot(10:10:50000,result12(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result22(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result32(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result42(1:5000,1),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('个流吞吐 (Gbps)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 40])
set(gca,'FontName','Times New Roman','FontSize',14);

%S0队列
subplot(6,4,6)
plot(10:10:50000,result12(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result22(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result32(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result42(1:5000,2),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S0队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 1000])
set(gca,'FontName','Times New Roman','FontSize',14);

%S1队列
subplot(6,4,7)
plot(10:10:50000,result12(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result22(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result32(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result42(1:5000,3),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S1队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 1000])
set(gca,'FontName','Times New Roman','FontSize',14);

%丢包
subplot(6,4,8)
plot(10:10:50000,result12(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result22(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result32(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result42(1:5000,4),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('累积丢包 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 max(100,5000)])
set(gca,'FontName','Times New Roman','FontSize',14);

%第三行
%单流吞吐
subplot(6,4,9)
plot(10:10:50000,result13(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result23(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result33(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result43(1:5000,1),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('个流吞吐 (Gbps)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 12.5])
set(gca,'FontName','Times New Roman','FontSize',14);

%S0队列
subplot(6,4,10)
plot(10:10:50000,result13(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result23(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result33(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result43(1:5000,2),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S0队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 500])
set(gca,'FontName','Times New Roman','FontSize',14);

%S1队列
subplot(6,4,11)
plot(10:10:50000,result13(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result23(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result33(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result43(1:5000,3),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S1队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 500])
set(gca,'FontName','Times New Roman','FontSize',14);

%丢包
subplot(6,4,12)
plot(10:10:50000,result13(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result23(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result33(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result43(1:5000,4),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('累积丢包 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 5])
set(gca,'FontName','Times New Roman','FontSize',14);


%第四行
%单流吞吐
subplot(6,4,13)
plot(10:10:50000,result14(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result24(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result34(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result44(1:5000,1),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('个流吞吐 (Gbps)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 40])
set(gca,'FontName','Times New Roman','FontSize',14);

%S0队列
subplot(6,4,14)
plot(10:10:50000,result14(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result24(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result34(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result44(1:5000,2),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S0队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 1000])
set(gca,'FontName','Times New Roman','FontSize',14);

%S1队列
subplot(6,4,15)
plot(10:10:50000,result14(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result24(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result34(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result44(1:5000,3),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S1队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 1000])
set(gca,'FontName','Times New Roman','FontSize',14);

%丢包
subplot(6,4,16)
plot(10:10:50000,result14(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result24(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result34(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result44(1:5000,4),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('累积丢包 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 max(100,50000)])
set(gca,'FontName','Times New Roman','FontSize',14);


%第五行
%单流吞吐
subplot(6,4,17)
plot(10:10:50000,result15(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result25(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result35(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result45(1:5000,1),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('个流吞吐 (Gbps)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 10])
set(gca,'FontName','Times New Roman','FontSize',14);

%S0队列
subplot(6,4,18)
plot(10:10:50000,result15(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result25(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result35(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result45(1:5000,2),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S0队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 500])
set(gca,'FontName','Times New Roman','FontSize',14);

%S1队列
subplot(6,4,19)
plot(10:10:50000,result15(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result25(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result35(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result45(1:5000,3),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S1队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 500])
set(gca,'FontName','Times New Roman','FontSize',14);

%丢包
subplot(6,4,20)
plot(10:10:50000,result15(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result25(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result35(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result45(1:5000,4),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('累积丢包 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 5])
set(gca,'FontName','Times New Roman','FontSize',14);


%第六行
%单流吞吐
subplot(6,4,21)
plot(10:10:50000,result16(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result26(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result36(1:5000,1),'LineWidth',2);
hold on;
plot(10:10:50000,result46(1:5000,1),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('个流吞吐 (Gbps)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 40])
set(gca,'FontName','Times New Roman','FontSize',14);

%S0队列
subplot(6,4,22)
plot(10:10:50000,result16(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result26(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result36(1:5000,2),'LineWidth',2);
hold on;
plot(10:10:50000,result46(1:5000,2),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S0队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 1000])
set(gca,'FontName','Times New Roman','FontSize',14);

%S1队列
subplot(6,4,23)
plot(10:10:50000,result16(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result26(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result36(1:5000,3),'LineWidth',2);
hold on;
plot(10:10:50000,result46(1:5000,3),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('S1队列 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 1000])
set(gca,'FontName','Times New Roman','FontSize',14);

%丢包
subplot(6,4,24)
plot(10:10:50000,result16(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result26(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result36(1:5000,4),'LineWidth',2);
hold on;
plot(10:10:50000,result46(1:5000,4),'LineWidth',2);
hold on;
xlabel('时间 (us)','FontName','Times New Roman','FontSize',14);
ylabel('累积丢包 (包)','FontName','Times New Roman','FontSize',14);
axis([0 50000 0 max(100,800000)])
set(gca,'FontName','Times New Roman','FontSize',14);
