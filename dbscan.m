%% DBSCAN
function [class,number]=dbscan(data)

% 定义参数Eps和MinPts
MinPts =2;
Eps = epsilon(data, MinPts);

[m,n] = size(data);%得到数据的大小

x = [(1:m)' data];
[m,n] = size(x);%重新计算数据集的大小
types = zeros(1,m);%用于区分核心点1，边界点0和噪音点-1
dealed = zeros(m,1);%用于判断该点是否处理过,0表示未处理过
dis = calDistance(x(:,2:n));
number = 1;%用于标记类

%% 对每一个点进行处理
for i = 1:m
    %找到未处理的点
    if dealed(i) == 0
        xTemp = x(i,:);
        D = dis(i,:);%取得第i个点到其他所有点的距离
        ind = find(D<=Eps);%找到半径Eps内的所有点
        
        %% 区分点的类型
        
        %边界点
        if length(ind) > 1 && length(ind) < MinPts+1
            types(i) = 0;   
            class(i) = 0;
        end
        %噪音点
        if length(ind) == 1
            types(i) = -1;
            class(i) = -1;
            dealed(i) = 1;
        end
        %核心点(此处是关键步骤)
        if length(ind) >= MinPts+1
            types(xTemp(1,1)) = 1;
            class(ind) = number;
            
            % 判断核心点是否密度可达
            while ~isempty(ind)   %判断ind是否不为空，不为空，则为1
                yTemp = x(ind(1),:);
                dealed(ind(1)) = 1;
                ind(1) = [];
                D = dis(yTemp(1,1),:);%找到与ind(1)之间的距离
                ind_1 = find(D<=Eps);
                
                if length(ind_1)>1 %处理非噪音点
                    class(ind_1) = number;
                    if length(ind_1) >= MinPts+1
                        types(yTemp(1,1)) = 1;
                    else
                        types(yTemp(1,1)) = 0;
                    end
                    
                    for j=1:length(ind_1)
                       if dealed(ind_1(j)) == 0
                          dealed(ind_1(j)) = 1;
                          ind=[ind ind_1(j)];   
                          class(ind_1(j))=number;
                       end                    
                   end
                end
            end
            number = number + 1;
        end
    end
end

% 最后处理所有未分类的点为噪音点
ind_2 = find(class==0);
class(ind_2) = -1;
types(ind_2) = -1;

% class
% number