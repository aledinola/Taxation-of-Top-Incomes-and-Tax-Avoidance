function A = load_NoAvoidanceExp_tran(dirname)
% PURPOSE
%load_NoAvoidanceExp_tran loads transitional dynamics result from folder 'dirname' and
%stores them into the structure A
% INPUTS
% dirname : name of the txt directory where the results of comparative statics 
%            are saved
% OUTPUTS
% A        : structure where variables are stored

if ~ischar(dirname)
    error('Input argument dirname must be a character')
end

%result = load(dirname);

% Read the number of periods on the transition path
tmp = load(fullfile(dirname,'KN_path.txt'));
T = numel(tmp);

% Transition path variables 
A.r = loadArray(fullfile(dirname,'r_path.txt'),[T,1]);
A.w = loadArray(fullfile(dirname,'wage_path.txt'),[T,1]);
A.Y = loadArray(fullfile(dirname,'Y_path.txt'),[T,1]);
A.K = loadArray(fullfile(dirname,'A_path.txt'),[T,1]); % A (asset) is equal to agg capital
A.avek = loadArray(fullfile(dirname,'avek_path.txt'),[T,4]); % first col: all entre
A.avek_allentre = A.avek(:,1);
A.Y_entre = loadArray(fullfile(dirname,'Y_entre_path.txt'),[T,1]);
A.taxes_inc = loadArray(fullfile(dirname,'taxes_inc_path.txt'),[T,1]);
A.taxes_corp = loadArray(fullfile(dirname,'taxes_corp_path.txt'),[T,1]);
A.taxes_div = loadArray(fullfile(dirname,'taxes_div_path.txt'),[T,1]);
A.taxes_ss = loadArray(fullfile(dirname,'taxes_ss_path.txt'),[T,1]);
A.tau_p = loadArray(fullfile(dirname,'tau_p_path.txt'),[T,1]);
A.tau_p = loadArray(fullfile(dirname,'tau_p_path.txt'),[T,1]);
A.share_entre = loadArray(fullfile(dirname,'share_entre_path.txt'),[T,4]); % first col: share of entre



end %end function "load_NoAvoidanceExp_tran"