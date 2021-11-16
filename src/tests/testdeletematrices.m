index2delete = find(in==0); %list R
TListLen = length(TList); %testl R
index2deleteLen = length(index2delete); %listl R
TemTList = zeros(4,4,TListLen-index2deleteLen);%output
TemqList = zeros(6,TListLen-index2deleteLen);
TemqListNormalized = zeros(6,TListLen-index2deleteLen);
TemxList = zeros(3,TListLen-index2deleteLen);
jumper = index2deleteLen;

for i = TListLen: -1: 1
    if (jumper > 0)
        if(~(i == index2delete(jumper)))
            TemTList(:,:,i-jumper) = TList(:,:,i);
            TemqList(:,i-jumper) = qList(:,i);
            TemqListNormalized(:,i-jumper) = qListNormalized(:,i);
            TemxList(:,i-jumper) = xList(:,i);
        else
            jumper = jumper - 1;
        end
    else
        TemTList(:,:,i) = TList(:,:,i);
        TemqList(:,i) = qList(:,i);
        TemqListNormalized(:,i) = qListNormalized(:,i);
        TemxList(:,i) = xList(:,i);
    end
end

TList = TemTList;
qList = TemqList;
qListNormalized = TemqListNormalized;
xList = TemxList;