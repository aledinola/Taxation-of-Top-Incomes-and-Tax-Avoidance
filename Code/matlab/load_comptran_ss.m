function A = load_comptran_ss(dirname)
% PURPOSE
%load_compstat loads steady state result from file 'filename' and
%stores them into the structure A
% INPUTS
% filename : name of the txt directory where the results of comparative statics 
%            are saved
% OUTPUTS
% A        : structure where variables are stored

if ~ischar(dirname)
    error('Input argument dirname must be a character')
end

A.Y = importdata(fullfile(dirname,'agg_Y.txt'));
A.C = importdata(fullfile(dirname,'agg_C.txt'));
A.K = importdata(fullfile(dirname,'agg_A.txt'));
A.r = importdata(fullfile(dirname,'r.txt'));
A.w = importdata(fullfile(dirname,'w.txt')); % this is not yet available
A.KN = importdata(fullfile(dirname,'KN.txt')); % this is not yet available
A.hsv_0 = importdata(fullfile(dirname,'hsv_0.txt'));
A.tau_p = importdata(fullfile(dirname,'tau_p.txt'));
% entrep. and LFO shares
A.share_entre = importdata(fullfile(dirname,'share_entre_act.txt'));
A.share_ep = importdata(fullfile(dirname,'share_EP_entre.txt'));
A.share_es = importdata(fullfile(dirname,'share_ES_entre.txt'));
A.share_ec = importdata(fullfile(dirname,'share_EC_entre.txt'));

if isfile(fullfile(dirname,'agg_Y_entre.txt'))
    A.Y_entre = importdata(fullfile(dirname,'agg_Y_entre.txt'));
end
if isfile(fullfile(dirname,'avek_all.txt'))
    A.agg_avek_allentre = importdata(fullfile(dirname,'avek_all.txt'));
end
if isfile(fullfile(dirname,'taxes_inc.txt'))
    A.taxes_inc = importdata(fullfile(dirname,'taxes_inc.txt'));
end
if isfile(fullfile(dirname,'taxes_corp.txt'))
    A.taxes_corp = importdata(fullfile(dirname,'taxes_corp.txt'));
end
if isfile(fullfile(dirname,'taxes_div.txt'))
    A.taxes_div = importdata(fullfile(dirname,'taxes_div.txt'));
end
if isfile(fullfile(dirname,'taxes_ss.txt'))
    A.taxes_ss = importdata(fullfile(dirname,'taxes_ss.txt'));
end



end