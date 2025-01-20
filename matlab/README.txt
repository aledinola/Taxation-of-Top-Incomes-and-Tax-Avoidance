
Matlab programs

matlab_plotter.m analyzes results from the benchmark steady state and plots policy functions etc.

NoAvoidanceExp_plotter.m analyzes results from no-avoidance experiments

compstat_plotter.m analyzes results from comparative statics of varying tau_h in steady-states. 
	Functions called:
	- load_compstat.m loads results from "compstat..." folders in "output" 
	- fun_plot_compstat.m make plots of various variables against tau_h
	
	
comptran_plotter.m analyzes results from varying tau_h with transition.
	Functions called:
	- load_comptran.m loads results from "comptran..." folders in "output"
	- fun_plot_comptran.m make plots, including CEV by tau_h, CEV by group given the optimal tau_h, and transition paths.
	
------------------
DICTIONARY
------------------

Figures produced by matlab_plotter are stored in subfolder "ss_bench/"
Fig name (ss_bench/...) Description											Variable Name
occpol_prob_W/EP/ES/EC	Occupation and legal form choice prob.				occpol_prob = muy./sum(muy,4)
phipol_ES/EC			Share of income declared as wage					phipol


Figures produced by NoAvoidanceExp_plotter are stored in subfolder "exp/"
Fig name (exp/...) 		Description											Variable Name
exp1/2/3/4_cev_qo		CEV by occupation and wealth bins					cev_qo
cev_z					CEV by occupation and LFO							cev_z


Figures produced by compstat_plotter (no fiscal neutraliy) have prefix "compstat_..."
Fig name (compstat_...) Description											Variable Name
taxes_inc_corp_div		Total tax revenue (inc., corp., and div. taxes)		taxes_inc_corp_div
base_taxes				Tax revenue by type, benchmark						taxes_inc and taxes_corp_div (base only, normalized)
share_entre_act			Share of entrepreneurs in young population			share_entre_act
share_lfo_top5_EP/S/C	Sole-prop./S/C share of entre., top 5% income		share_lfo_top5(:,1/2/3)
share_EP/S/C_entre		Sole-proprietor/S/C share of entrepreneurs			share_EP/S/C_entre
base_share_wage			Share of income declared as wage, benchmark			share_wage_ES/C
Y/C/K/N					Total output/consumption/capital/labor				Y/C/K/N
r/w						Interest rate/wage									r/w
	
Figures produced by comptran_plotter (with fiscal neutraliy) have prefix "comptran_..."
Fig name (comptran_...) Description											Variable Name
cev_vec					Aggregate CEV										cev_vec
cev_aggcomp_vec			Aggregate component of CEV							cev_aggcomp_vec
cev_distcomp_vec		Distributional component of CEV						cev_distcomp_vec
base/CF2_cev_qo			CEV by group, optimal tau_h, benchmark/no-avoidance	cev_qo
cev_z					CEV by occupation and LFO, optimal tau_h			cev_z
Y/C/N/K/r/w_path		Transition Paths: total output/consumption/labor 	Y/C/N/K/r/w_path
							supply/capital/interest rate/wage
base/CF2_shares_path	Transition paths: occ. shares						share_entre/ep/es/ec, normalized
											

