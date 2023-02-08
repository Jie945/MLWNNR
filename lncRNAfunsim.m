function [mms]=lncRNAfunsim(DDS,M)
mrows=size(M,1);
mms=zeros(mrows,mrows);
for i=1:mrows
    index1=find(M(i,:)==1);
    if length(index1)==0
        mms(i,:)=0;
    else
        for j=1:mrows
            index2=find(M(j,:)==1);
            if length(index2)==0
                mms(i,j)=0;
            else
                sim1=zeros(1,length(index1));
                for m=1:length(index1)
                    sim1(m)=max(DDS(index2,index1(m))');
                end
                sim2=zeros(1,length(index2));
                for n=1:length(index2)
                    sim2(n)=max(DDS(index1,index2(n))');
                end
                mms(i,j)=(sum(sim1)+sum(sim2))/(length(index1)+length(index2));
            end
        end
    end
end
for k=1:mrows
    mms(k,k)=1;
end
end




