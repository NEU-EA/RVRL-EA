function [state,rho] = getState(PopObj,V)
% ��õ�ǰ������״̬

    PopObj        = PopObj - repmat(min(PopObj,[],1),size(PopObj,1),1);
    Front = NDSort(PopObj,size(PopObj,1));
    [rho,pi]   = Associate(PopObj,V); %Pi��ʾÿ�����������������˭��rho��ʾÿ�����������м�������   
    state = ones(1,size(V,1));
    for i = 1:size(V,1)
        if rho(i)==0
            state(i) = 2;%��������û�н�
        end
        if rho(i) == 1
            if Front(pi==i)==1
                state(i) = 1;%һ����֧���
            else
                state(i) = 4;%һ��֧���
            end
        end
        if rho(i) >1
            if Front(pi==i)==1
                state(i) = 3;%���������ж����֧���
            else
                if Front(pi==i)>1
                     state(i) = 5;%���������ж��֧���
                else
                    state(i) = 6;%���������ж����Ͻ�
                end
            end
        end
    end
%     state(rho==1) = 1;
%     state(rho==0) = 2;
%     state(rho>1) = 3;
%     inValid       = find(rho==0);
%     R = getReward();
%     q=ReinforcementLearningUpdateR(R);
%     V(inValid,:)  = rand(length(inValid),size(V,2)).*repmat(max(PopObj,[],1),length(inValid),1);
end

function [rho,pi] = Associate(PopObj,Z)
% Associate each solution with one reference point

    %% Calculate the distance of each solution to each reference vector
    NormP    = sqrt(sum(PopObj.^2,2));
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(NormP,1,size(Z,1)).*sqrt(1-Cosine.^2);
    
    %% Associate each solution with its nearest reference point
    [~,pi] = min(Distance',[],1);
    
    %% Calculate the number of associated solutions of each reference point
    rho = hist(pi,1:size(Z,1));    
end
