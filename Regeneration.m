function V_new = Regeneration(V_old,y) % V_old��ʾ�ο�������y��ʾ�����Ŀ��ֵ
%REGENERATION ����EoDģ���������ɲο�����
%   1���ֲ�ģ�ͣ�ͨ�����࣬���������࣬Ȼ��ÿ���๹��һ��ģ��
%   2��ȫ��ģ�ͣ�����ȫ�������ݹ���һ��ģ��
%   3��ÿ���Ӿֲ���ȫ�����е�ģ��������ѡ��һ��ģ�Ͳ���
[~,rho]=getState(y,V_old);
V_new=zeros(1,size(V_old,2));

%% Step 1.Clustering
V_old = V_old(rho~=0,:);
[class,~]=dbscan(V_old);% �������ÿһ�д���һ���ο�����
% �������classΪm*1�����飬��ʾm���ο��������ڵľ������
[~,n]=size(V_old);% n��ʾĿ����
N_W= max(class); %length(unique(class));% Number of clusters
if all(class == -1)
    N_W=0;
end
    
% ��һ��y
% y_norm=zeros(size(y));% ���������
% for j=1:n
%     y_norm(:,j)=mapminmax(y(:,j)',0,1)';
% end

%% EoD learning model construction
for i=1:n % Ŀ������n
    clear PRO;% �������PRO
    clear r_matrix;% �������r_matrix
    % local probability model
    for k=1:N_W % ��������N_W
        V_ind=find(class==k);% ��i�������������������вο����������
        
        W=length(V_ind);% �����вο�����������
        
        [C,r_arr]=Cij(i,W,V_old(V_ind,:));
        r_matrix{k}=r_arr;
        % ��ֲ�����
        PRO{k}=C./sum(C,2);% Ԫ�����飬�洢��ͬ����ĸ�������
        
    end
    % global probability model
    W=size(V_old,1);
    [C,r_arr]=Cij(i,W,V_old);
    r_matrix{N_W+1}=r_arr;
    % ��ȫ�ָ���
    PRO{N_W+1}=C./sum(C,2);% Ԫ�����飬�洢��ͬ����ĸ�������
    
    %% Sampling
    count=size(PRO,2);% ����count������ģ��
    index=floor(rand*count)+1;% ���ѡȡģ��
    
    % ���̶�ȷ���²ο������ĵ�iά��ֵ
    pro=PRO{index};
    r_temp=r_matrix{index};
    pro=pro*triu(ones(size(pro,2),size(pro,2)),0);% ���ۻ�����
    rand_temp=rand;
    for z=1:size(pro,2)
        if rand_temp<pro(z)
            V_temp=rand*(r_temp(z+1)-r_temp(z))+r_temp(z);
            break;
        end
    end
    V_new(1,i)=V_temp;
end
% PRO_i PRO_0 �ܹ�N_W+1��ֱ��ͼ����ģ�ͣ������Թ�ʽ(11)����ʽ�����飬ÿ��ģ�ͷ�ΪW+1��bin;
% W�ǲο������ĸ���(�����Ǿ���k�еĲο�����������Ҳ������ȫ�ֵĲο���������).
% r����һ�б�ʾһ�������bin���ֻ���һ��ȫ��bin���֣���λ�洢ri_��ri-
% ͨ�����̶ĵķ�ʽѡȡ��j��bin�����µĲο�����w'
end

%% ��ÿ��bin�ĸ߶�
function [C,r_arr]=Cij(i,W,y_norm)
% ȷ��ri_��ri-
r_min=0;
r_max=1;

% ���ֳ�W��bin
r_arr=linspace(r_min,r_max,W+1); % ������

% ��Ci��j
C=zeros(1,W);

for q=1:size(y_norm,1)% ��Ч����������
    for h=2:W+1 % ���ֵ�bin����
        if y_norm(q,i)<r_arr(1,h) %&& y_norm(q,i)>r_arr(1,h-1) 
            C(1,h-1)=C(1,h-1)+1;% ͳ�Ƹ�������h-1bin�еĸ���
            break;
        end
    end
end

C(1,find(C==0))=0.1; % �߶�Ϊ0��bin��һ����С��ֵ
end