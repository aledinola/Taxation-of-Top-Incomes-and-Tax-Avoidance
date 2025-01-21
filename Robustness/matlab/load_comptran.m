function A = load_comptran(dirname)
% PURPOSE
%load_compstat loads comparative statics result from file 'filename' and
%stores them into the structure A
% INPUTS
% filename : name of the txt directory where the results of comparative statics 
%            are saved
% OUTPUTS
% A        : structure where variables are stored

if ~ischar(dirname)
    error('Input argument dirname must be a character')
end

%result = load(dirname);

% Read the size of the tau_h vector
tmp = load(fullfile(dirname,'tau_h_vec.txt'));
n_tau_h = numel(tmp);

% Read the number of periods on the transition path
tmp = load(fullfile(dirname,'KN_path.txt'));
T = numel(tmp);

% Read the number of percentile groups for group CEV analysis
tmp = load(fullfile(dirname,'cev_q.txt'));
nq = numel(tmp);

%% Fields 
% Vectors of tau_h and CEV
A.tau_h        = loadArray(fullfile(dirname,'tau_h_vec.txt'),[n_tau_h,1]);
A.cev_vec       = loadArray(fullfile(dirname,'cev_vec.txt'),[n_tau_h,1]);
if isfile(fullfile(dirname,'cev_aggcomp_vec.txt'))
    A.cev_aggcomp_vec = loadArray(fullfile(dirname,'cev_aggcomp_vec.txt'),[n_tau_h,1]);
end
if isfile(fullfile(dirname,'cev_distcomp_vec.txt'))
    A.cev_distcomp_vec = loadArray(fullfile(dirname,'cev_distcomp_vec.txt'),[n_tau_h,1]);
end

% Transition path variables (transition to the optimal tax rate)
A.KN = loadArray(fullfile(dirname,'KN_path.txt'),[T,1]);
A.hsv_0 = loadArray(fullfile(dirname,'hsv_0_path.txt'),[T,1]);
A.tau_p = loadArray(fullfile(dirname,'tau_p_path.txt'),[T,1]);
A.Y = loadArray(fullfile(dirname,'Y_path.txt'),[T,1]);
A.K = loadArray(fullfile(dirname,'A_path.txt'),[T,1]); % A (asset) is equal to agg capital
A.C = loadArray(fullfile(dirname,'C_path.txt'),[T,1]);
A.r = loadArray(fullfile(dirname,'r_path.txt'),[T,1]);
A.w = loadArray(fullfile(dirname,'wage_path.txt'),[T,1]);
A.share = loadArray(fullfile(dirname,'share_entre_path.txt'),[T,4]); % columns: all entre, EP, ES, EC
A.share(T,:) = A.share(T-1,:); % last period is undefined.

A.share_entre = A.share(:,1);
A.share_ep = A.share(:,2);
A.share_es = A.share(:,3);
A.share_ec = A.share(:,4);


% Other CEV measures at the optimal tau_h
A.cev_q = loadArray(fullfile(dirname,'cev_q.txt'),[nq,1]);
A.cev_qo = loadArray(fullfile(dirname,'cev_qo.txt'),[nq,2]); % first column is worker, second column is entre
A.cev_z = loadArray(fullfile(dirname,'cev_z.txt'),[4,1]); % Worker/EP/ES/EC

end %end function "load_comptran"