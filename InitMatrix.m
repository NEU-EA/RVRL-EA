% 初始化Q-learning的q-table

function [q,q1,R]=InitMatrix(N)

%% RLEA第一版
%     R=cell(1,N);
%     q=cell(1,N);
%     q1=cell(1,N);
%     for i=1:N
%         R{i} = ones(3,3);
%         q{i} = zeros(3,3);
%         q1{i} = ones(3,3)*inf;
%     end
    
%% RLEA第二版
% R = [-inf,0;
%     0 0;
%     0 0;
%     0 -inf];
% q=zeros(size(R));        % initialize Q as zero
% q1=ones(size(R))*inf;    % initialize previous Q as big number

%% RLEA终版
    R=cell(1,N);
    q=cell(1,N);
    q1=cell(1,N);
   for i=1:N
        R{i} = ones(6,4);
        q{i} = zeros(6,4);
        q1{i} = ones(6,4)*inf;
    end
end
