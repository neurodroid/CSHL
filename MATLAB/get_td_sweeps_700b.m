%puts data from 700b with inputs in AI0 and AI1 into a single matrix where
%each column in response and cmd are an individual sweep

function [response,cmd,si] = get_td_sweeps_700b(filename) 


my_tdms_struct = TDMS_getStruct(filename);

Ft = str2num(my_tdms_struct.Props.Sampling_Rate_AI_); %sample rate
si = 1000/str2double(my_tdms_struct.Props.Sampling_Rate_AI_);%sample interval ms


%%coalate data

fn = fieldnames(my_tdms_struct.AI0);
sweeps = length(fn)-2;%first two field names are descriptors, rest ar sweeps
for i = 1:sweeps
    response(:,i) = eval(strcat('my_tdms_struct.AI0.',cell2mat(fn(i+2)),'.data'));
    cmd(:,i) = eval(strcat('my_tdms_struct.AI1.',cell2mat(fn(i+2)),'.data'));

end

end