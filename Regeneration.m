function V_new = Regeneration(V_old,y) % V_old表示参考向量，y表示个体的目标值
%REGENERATION 根据EoD模型重新生成参考向量
%   1，局部模型，通过聚类，把向量分类，然后每个类构造一个模型
%   2，全局模型，即用全部的数据构造一个模型
%   3，每个从局部和全局所有的模型中任意选择一个模型采样
[~,rho]=getState(y,V_old);
V_new=zeros(1,size(V_old,2));

%% Step 1.Clustering
V_old = V_old(rho~=0,:);
[class,~]=dbscan(V_old);% 输入变量每一行代表一个参考向量
% 输出变量class为m*1的数组，表示m个参考向量所在的聚类序号
[~,n]=size(V_old);% n表示目标数
N_W= max(class); %length(unique(class));% Number of clusters
if all(class == -1)
    N_W=0;
end
    
% 归一化y
% y_norm=zeros(size(y));% 生成零矩阵
% for j=1:n
%     y_norm(:,j)=mapminmax(y(:,j)',0,1)';
% end

%% EoD learning model construction
for i=1:n % 目标数：n
    clear PRO;% 清除变量PRO
    clear r_matrix;% 清除变量r_matrix
    % local probability model
    for k=1:N_W % 聚类数：N_W
        V_ind=find(class==k);% 第i个聚类中所包含的所有参考向量的序号
        
        W=length(V_ind);% 聚类中参考向量的数量
        
        [C,r_arr]=Cij(i,W,V_old(V_ind,:));
        r_matrix{k}=r_arr;
        % 求局部概率
        PRO{k}=C./sum(C,2);% 元胞数组，存储不同聚类的概率数组
        
    end
    % global probability model
    W=size(V_old,1);
    [C,r_arr]=Cij(i,W,V_old);
    r_matrix{N_W+1}=r_arr;
    % 求全局概率
    PRO{N_W+1}=C./sum(C,2);% 元胞数组，存储不同聚类的概率数组
    
    %% Sampling
    count=size(PRO,2);% 共有count个概率模型
    index=floor(rand*count)+1;% 随机选取模型
    
    % 轮盘赌确定新参考向量的第i维的值
    pro=PRO{index};
    r_temp=r_matrix{index};
    pro=pro*triu(ones(size(pro,2),size(pro,2)),0);% 求累积概率
    rand_temp=rand;
    for z=1:size(pro,2)
        if rand_temp<pro(z)
            V_temp=rand*(r_temp(z+1)-r_temp(z))+r_temp(z);
            break;
        end
    end
    V_new(1,i)=V_temp;
end
% PRO_i PRO_0 总共N_W+1个直方图概率模型，都是以公式(11)的形式的数组，每个模型分为W+1个bin;
% W是参考向量的个数(可以是聚类k中的参考向量个数，也可以是全局的参考向量个数).
% r矩阵一行表示一个聚类的bin划分或者一个全局bin划分，首位存储ri_和ri-
% 通过轮盘赌的方式选取第j个bin生成新的参考向量w'
end

%% 求每个bin的高度
function [C,r_arr]=Cij(i,W,y_norm)
% 确定ri_和ri-
r_min=0;
r_max=1;

% 划分成W个bin
r_arr=linspace(r_min,r_max,W+1); % 行向量

% 求Ci，j
C=zeros(1,W);

for q=1:size(y_norm,1)% 有效向量的数量
    for h=2:W+1 % 划分的bin数量
        if y_norm(q,i)<r_arr(1,h) %&& y_norm(q,i)>r_arr(1,h-1) 
            C(1,h-1)=C(1,h-1)+1;% 统计个体落在h-1bin中的个数
            break;
        end
    end
end

C(1,find(C==0))=0.1; % 高度为0的bin赋一个很小的值
end