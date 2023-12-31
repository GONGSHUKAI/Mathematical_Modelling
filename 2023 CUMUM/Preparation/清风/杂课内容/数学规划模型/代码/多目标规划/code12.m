%%  多目标规划问题
w1 = 0.4;  w2 = 0.6;  % 两个目标函数的权重  x1 = 5  x2 = 2
w1 = 0.5;  w2 = 0.5;  % 两个目标函数的权重  x1 = 5  x2 = 2
w1 = 0.3;  w2 = 0.7;  % 两个目标函数的权重  x1 = 1  x2 = 6
c = [w1/30*2+w2/2*0.4 ;w1/30*5+w2/2*0.3];  % 线性规划目标函数的系数
A = [-1 -1];  b = -7; % 不等式约束
lb = [0 0]'; ub = [5 6]'; % 上下界
[x,fval] = linprog(c,A,b,[],[],lb,ub)
f1 = 2*x(1)+5*x(2)
f2 = 0.4*x(1) + 0.3*x(2)


%% 敏感性分析
clear;clc
W1 = 0.1:0.001:0.5;  W2 = 1- W1;  
n =length(W1);
F1 = zeros(n,1);  F2 = zeros(n,1);   X1 = zeros(n,1);  X2 = zeros(n,1);   FVAL = zeros(n,1);
A = [-1 -1];  b = -7; % 不等式约束
lb = [0 0]; ub = [5 6]; % 上下界
for i = 1:n
    w1 = W1(i);  w2 = W2(i);
    c = [w1/30*2+w2/2*0.4 ;w1/30*5+w2/2*0.3];  % 线性规划目标函数的系数
    [x,fval] = linprog(c,A,b,[],[],lb,ub);
    F1(i) = 2*x(1)+5*x(2);
    F2(i) = 0.4*x(1) + 0.3*x(2);
    X1(i) = x(1);
    X2(i) = x(2);
    FVAL(i) = fval;
end

% 「Matlab」“LaTex字符汇总”讲解：https://blog.csdn.net/Robot_Starscream/article/details/89386748
% 在图上可以加上数据游标，按住Alt加鼠标左键可以设置多个数据游标出来。
figure(1) 
plot(W1,F1,W1,F2)
xlabel('f_{1}的权重') 
ylabel('f_{1}和f_{2}的取值')
legend('f_{1}','f_{2}')

figure(2)
plot(W1,X1,W1,X2)
xlabel('f_{1}的权重') 
ylabel('x_{1}和x_{2}的取值')
legend('x_{1}','x_{2}')

figure(3)
plot(W1,FVAL)  % 看起来是两个直线组合起来的下半部分
xlabel('f_{1}的权重') 
ylabel('综合指标的值')



% % 注意：代码文件仅供参考，一定不要直接用于自己的数模论文中
% % 国赛对于论文的查重要求非常严格，代码雷同也算作抄袭
% % 视频中提到的附件可在售后群（购买后收到的那个无忧自动发货的短信中有加入方式）的群文件中下载。包括讲义、代码、我视频中推荐的资料等。
% % 关注我的微信公众号《数学建模学习交流》，后台发送“软件”两个字，可获得常见的建模软件下载方法；发送“数据”两个字，可获得建模数据的获取方法；发送“画图”两个字，可获得数学建模中常见的画图方法。另外，也可以看看公众号的历史文章，里面发布的都是对大家有帮助的技巧。
% % 购买更多优质精选的数学建模资料，可关注我的微信公众号《数学建模学习交流》，在后台发送“买”这个字即可进入店铺(我的微店地址：https://weidian.com/?userid=1372657210)进行购买。
% % 视频价格不贵，但价值很高。单人购买观看只需要58元，三人购买人均仅需46元，视频本身也是下载到本地观看的，所以请大家不要侵犯知识产权，对视频或者资料进行二次销售。
% % 如何修改代码避免查重的方法：https://www.bilibili.com/video/av59423231（必看）
