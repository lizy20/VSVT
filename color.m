
li=['#009392'; '#39B185'; '#9CCB86'; '#E9E29C'; '#EEB479'; '#E88471'; '#CF597E'];
for i =1:7 
    num=sixteen2ten(li(i,:))
end

function num =sixteen2ten(string)
    exchange_list='0123456789ABCDEF#';
    num=zeros(1,3);
    for i=1:3
        tempCoe1=find(exchange_list==string(i*2))-1;
        tempCoe2=find(exchange_list==string(i*2+1))-1;
        num(i)=16*tempCoe1+tempCoe2;
    end
end

function string =ten2sixteen(num)
%the num should be a 1x3 Integer mat limited in [0 255]
exchange_list={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
string='#';
for i=1:3
    temp_num=num(i);
    string(1+i*2-1)=exchange_list{(temp_num-mod(temp_num,16))/16+1};
    string(1+i*2)=exchange_list{mod(temp_num,16)+1};
end
end
% c=[ 0   147   146; 57   177   133;156   203   134;233   226   156;238   180   121;232   132   113;
   % 207    89   126];
% cmap=[ 057,081,162;  114,170,207; 202,232,242; 254,251,186;253,185,107;236,093,059;...
% 168,003,038]/225;
