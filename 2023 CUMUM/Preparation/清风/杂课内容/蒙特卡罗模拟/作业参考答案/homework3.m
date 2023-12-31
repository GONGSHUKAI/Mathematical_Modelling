%%  蒙特卡罗求解非线性规划问题
% min f(x) =2*(x1^2)+x2^2-x1*x2-8*x1-3*x2
% s.t.
% (1) 3*x1+x2>9
% (2) x1+2*x2<16
% (3) x1>0 & x2>0

%% （1）初次寻找最小值的代码
clc,clear;
format long g   %可以将Matlab的计算结果显示为一般的长数字格式（默认会保留四位小数，或使用科学计数法）
tic %计算tic和toc中间部分的代码的运行时间
n=10000000; %生成的随机数组数
x1=unifrnd(0,16,n,1);  % 生成在[0,16]之间均匀分布的随机数组成的n行1列的向量构成x1
x2=unifrnd(0,8,n,1);  % 生成在[0,8]之间均匀分布的随机数组成的n行1列的向量构成x2
fmin=+inf; % 初始化函数f的最小值为正无穷（后续只要找到一个比它小的我们就对其更新）
for i=1:n
    x = [x1(i), x2(i)];  %构造x向量, 这里千万别写成了：x =[x1, x2]
    if (3*x(1)+x(2)>9)  &  (x(1)+2*x(2)<16)     % 判断是否满足条件
        result = 2*(x(1)^2)+x(2)^2-x(1)*x(2)-8*x(1)-3*x(2);  % 如果满足条件就计算函数值
        if  result  < fmin  % 如果这个函数值小于我们之前计算出来的最小值
            fmin = result;  % 那么就更新这个函数值为新的最小值
            X = x;  % 并且将此时的x1 x2 保存到相应的变量中
        end
    end
end
disp(strcat('蒙特卡罗模拟得到的最小值为',num2str(fmin)))
disp('最小值处x1 x2的取值为：')
disp(X)
toc %计算tic和toc中间部分的代码的运行时间


%% （2）缩小范围重新模拟得到更加精确的取值
clc,clear;
tic %计算tic和toc中间部分的代码的运行时间
n=10000000; %生成的随机数组数
x1=unifrnd(2,3,n,1);  % 生成在[2,3]之间均匀分布的随机数组成的n行1列的向量构成x1
x2=unifrnd(2,3,n,1);  % 生成在[2,3]之间均匀分布的随机数组成的n行1列的向量构成x2
fmin=+inf; % 初始化函数f的最小值为正无穷（后续只要找到一个比它小的我们就对其更新）
for i=1:n
    x = [x1(i), x2(i)];  %构造x向量, 这里千万别写成了：x =[x1, x2]
    if (3*x(1)+x(2)>9)  &  (x(1)+2*x(2)<16)     % 判断是否满足条件
        result = 2*(x(1)^2)+x(2)^2-x(1)*x(2)-8*x(1)-3*x(2);  % 如果满足条件就计算函数值
        if  result  < fmin  % 如果这个函数值小于我们之前计算出来的最小值
            fmin = result;  % 那么就更新这个函数值为新的最小值
            X = x;  % 并且将此时的x1 x2 保存到相应的变量中
        end
    end
end
disp(strcat('蒙特卡罗模拟得到的最小值为',num2str(fmin)))
disp('最小值处x1 x2的取值为：')
disp(X)
toc %计算tic和toc中间部分的代码的运行时间




% % 注意：代码文件仅供参考，一定不要直接用于自己的数模论文中
% % 国赛对于论文的查重要求非常严格，代码雷同也算作抄袭
% % 视频中提到的附件可在售后群（购买后收到的那个无忧自动发货的短信中有加入方式）的群文件中下载。包括讲义、代码、我视频中推荐的资料等。
% % 关注我的微信公众号《数学建模学习交流》，后台发送“软件”两个字，可获得常见的建模软件下载方法；发送“数据”两个字，可获得建模数据的获取方法；发送“画图”两个字，可获得数学建模中常见的画图方法。另外，也可以看看公众号的历史文章，里面发布的都是对大家有帮助的技巧。
% % 购买更多优质精选的数学建模资料，可关注我的微信公众号《数学建模学习交流》，在后台发送“买”这个字即可进入店铺(我的微店地址：https://weidian.com/?userid=1372657210)进行购买。
% % 视频价格不贵，但价值很高。单人购买观看只需要58元，三人购买人均仅需46元，视频本身也是下载到本地观看的，所以请大家不要侵犯知识产权，对视频或者资料进行二次销售。
% % 如何修改代码避免查重的方法：https://www.bilibili.com/video/av59423231（必看）