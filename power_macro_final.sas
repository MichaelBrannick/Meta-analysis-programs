
%Macro Metapower (test= , model= , raw_data= , alpha= , tau2= , heterogeneity = , 
  n1= , n2= , k = , eff_type = , T = , Dataset= , B= , v= , x= , es= , p= , weight= ) ;
* +-------------------------------------------------------------------------------------------------+
   macro arguments:
Test= Specifies the type of test for which power calculations are desired: M (mean), QT (overall heterogeneity), contrast (group contrasts), QB (categorical moderator), QW (heterogeneity across levels of a categorical variable), Reg (omnibus test and regression coefficients), or QE (heterogeneity in regression).
Model= Specifies the type of model, fixed or random. Random also refers to a mixed model.
Raw_Data= Specifies whether power analysis is based on raw data, yes, or not, no.   
Alpha = Specifies the nominal alpha level.
Tau2 = Specifies the amount of between study variance.  A value of 99 will allow you to specify a heterogeneity ratio instead of this value. 
Heterogeneity = This is the ratio of between to within variability discussed by Hedges and Pigott (2001,2004).
n1= / n2 = Specifies average sample size for the collection of studies. Both arguments are relevant for d, but only n1= is evaluated for the effect sizes z and r. For the odds ratio, specify the variance of the log odds in the average study.
Eff_type= Specify the type of effect size being used: standardized mean difference (d), correlation coefficient (r), Fishers z transformation of the correlation coefficient (z), or the log of the odds ratio (or).  
T= Specify the average magnitude of the effect size in the effect size units specified in the Eff_type= argument.   
Dataset= Specify a SAS dataset.
B= Specify a column from an external dataset that corresponds to standardized regression coefficients. 
V= Specify one or more columns from an external dataset that contain the variances of the study effect sizes.
X= Specify one or more columns from an external dataset that contain the values of predictor variables in regression. This corresponds to the design matrix. 
ES= Specify one or more columns from an external dataset that contain the effect sizes from each individual study.	
P= Corresponds to the number of categories in a test of residual variability in a categorical model or the number of predictors in a test of residual variability in regression.
Weight= Specify a column from an external dataset that corresponds to weights assigned to groups in a test of contrasts. 
  +-------------------------------------------------------------------------------------------------+;

proc iml;
* +-----------------------------------+
    Test of mean effect size (M)= 0
  +-----------------------------------+;
if &test = 'M' then do;
	if &raw_data= 'no' then do;
		if &eff_type = 'z' then do;
			z = .5 * log((1+&T) / (1-&T));
			common_var = 1 / (&n1 - 3);
			var_overall = common_var / &k;
		end;
		if &eff_type=  'r' then do;
			common_var=(((1-&T**2)**2)/(&n1-1));
			var_overall = common_var / &k;
		end; 
		if &eff_type = 'd' then do;
			common_var = ((&n1+&n2)/(&n1*&n2))+((&T**2)/(2*(&n1+&n2)));
			var_overall = common_var / &k;
		end;
		if &eff_type = 'or' then do;
			common_var=&n1;
			var_overall = common_var / &k;
		end;
	
		if &model = 'fixed' then do;
			power = 1-probnorm((probit (1-(&alpha/2)))-(&T/(sqrt (var_overall))))+ probnorm((-probit (1-(&alpha/2)))-(&T/(sqrt (var_overall))));
		end;

		if &model = 'random' then do;
			if &tau2 ^= 99 then do; 
				power = 1-probnorm((probit (1-(&alpha/2)))-(&T/(sqrt (((&tau2)+ common_var)/&k)))) + probnorm((-probit (1-(&alpha/2)))-(&T/(sqrt (((&tau2)+ common_var)/&k))));
			end;
			if &tau2 = 99 then do;
				power = 1-probnorm((probit (1-(&alpha/2)))-(&T/(sqrt (((&heterogeneity*common_var)+ common_var)/&k))))+ probnorm((-probit (1-(&alpha/2)))-(&T/(sqrt (((&heterogeneity*common_var)+ common_var)/&k))));
			end;
		end;
	end;
	if &raw_data= 'yes' then do;
		use &Dataset;
		read all var{&es} into es;
		read all var{&v}  into v; 
		k = nrow(es);
		df = k - 1 ;
		w= 1/v;
		var_overall=1/sum(w);
		
		if &model = 'fixed' then do;
			/*var_overall=1/sum(w);
			tau2=99;*/
			power = 1-probnorm((probit (1-(&alpha/2)))-(&T/(sqrt (var_overall))))+ probnorm((-probit (1-(&alpha/2)))-(&T/(sqrt (var_overall))));
		end;
		if &model = 'random' then do;
			mes = sum(es#w)/sum(w);
			q= sum(w#(es-mes)##2);
			tau2= (q - df)/(sum(w)-(sum(w#w)/sum(w)));
				if tau2<0 then do ; 
				tau2= 0 ;
				end ;
				if tau2>= 0 then do;
				tau2=tau2;
				end; 
			v_random = (v + tau2);
			w_random= 1/v_random;
			var_random_overall = 1/(sum(w_random));
			power =  1-probnorm((probit (1-(&alpha/2)))-(&T/(sqrt (var_random_overall))))+ probnorm((-probit (1-(&alpha/2)))-(&T/(sqrt (var_random_overall))));
		end;
	end;
end;
* +--------------------------------------------+
   Test of heterogeneity in effect sizes (QT)= 0
  +--------------------------------------------+;
if &test= 'QT' then do;
	if &raw_data= 'no' then do;
		if &eff_type = 'z' then do;
			z = .5 * log((1+&T) / (1-&T));
			common_var = 1 / (&n1 - 3);
			var_overall = common_var / &k;
		end;
		if &eff_type=  'r' then do;
			common_var=(((1-&T**2)**2)/(&n1-1));
			var_overall = common_var / &k;
		end; 
		if &eff_type = 'd' then do;
			common_var = ((&n1+&n2)/(&n1*&n2))+((&T**2)/(2*(&n1+&n2)));
			var_overall = common_var / &k;
		end;
		if &eff_type = 'or' then do;
			common_var=&n1;
			var_overall = common_var / &k;
		end;
		if &model= 'fixed' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (&k-1)), (&k-1), (&tau2/common_var)*(&k-1));
			end;
			if &tau2 = 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (&k-1)), (&k-1), &heterogeneity*(&k-1));
			end;
		end;		
		if &model= 'random' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (((cinv ((1-&alpha), (&k-1))*common_var)/(common_var+ &tau2)), (&k-1), 0);
			end;
			if &tau2 =99 then do;
				power= 1-probchi (((cinv ((1-&alpha), (&k-1))*common_var)/(common_var+ (common_var*&heterogeneity))), (&k-1), 0);
			end;
		end;
	end;
	if &raw_data= 'yes' then do;
		use &Dataset;
		read all var{&es} into es;
		read all var{&v}  into v; 
		k = nrow(es);
		df = k - 1 ;
		w= 1/v;
		var_overall=1/sum(w);
		average_var=var_overall*k;
		if &model= 'fixed' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (df)), df, (&tau2/average_var)*(df));
			end;
			if &tau2 = 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (df)), df, ((&heterogeneity*average_var)/average_var)*(df));
			end;
		end;
		/*if &model= 'random' then do;
			sum_w= sum(w);
			sum_ww=sum(w##2);
			sum_www=sum(w##3);
			a= sum_w-(sum_ww/sum_w);	
				if &tau2 ^= 99 then do;
					mean=(a*&tau2)+df; 
					variance= (2*df)+(4*a*&tau2)+(2*(sum_ww-(sum_www/sum_w)+((sum_ww##2)/(sum_w##2))))*(&tau2##2);
				end;
			if &tau2= 99 then do;
				mean=(a*(&heterogeneity*average_var))+(df); 
				variance= (2*df)+(4*a*(&heterogeneity*average_var))+(2*(sum_ww-(sum_www/sum_w)+((sum_ww##2)/(sum_w##2))))*((&heterogeneity*average_var)##2);
			end;
				r=(mean##2)/variance;
				m= mean/variance;
				power=1-probgam ((cinv ((1-&alpha), df))/r, m);*something wrong here;*/
		if &model= 'random' then do;
			sum_w= sum(w);
			sum_ww=sum(w##2);
			sum_www=sum(w##3);
			a= sum_w-(sum_ww/sum_w);
			if &tau2 ^= 99 then do;
				mean=(a*&tau2)+(df); 
				b=(df)+(2*a*&tau2)+(sum_ww-(2*(sum_www/sum_w))+((sum_ww##2)/(sum_w##2)))*(&tau2##2);
			end;
			if &tau2= 99 then do;
				mean=(a*(&heterogeneity*average_var))+(df); 
				b=(df)+(2*a*(&heterogeneity*average_var))+(sum_ww-(2*(sum_www/sum_w))+((sum_ww##2)/(sum_w##2)))*((&heterogeneity*average_var)##2);
			end;
				variance=2*b;
				r=variance/(2*mean);
				s=(2*(mean##2))/variance;
				power=1-probchi (cinv ((1-&alpha), df)/r, s, 0);	
		end;
	end;	
end;


* +------------------------------------------------------------+
   Test of fixed paramters in regression model (QR and Betas)= 0 
  +------------------------------------------------------------+;
if &test= 'Reg' then do;
 use &Dataset;
 read all var{&B} into B_vec;
 read all var{&x} into x;
 read all var{&v} into var_vec;
 read all var{&es}  into es;
 k= ncol(x);*p;
 n= nrow(x);*k;
 design = j(n,1,1) || x ;
 B_vec_remove= remove (B_vec,((k+1):n));
 
if &model= 'fixed' then do;
	W=1/var_vec;
	W=diag(W);
	B_cov= design`*W*design;
	B_cov=inv(B_cov);
	submatrix=B_cov[(k-1):(k+1),(k-1):(k+1)];
	B_SE=sqrt(vecdiag(submatrix));
	B_cov_inv= inv(submatrix);
	QR= B_vec_remove*B_cov_inv*B_vec_remove`;
	QR_power= 1-probchi (cinv ((1-&alpha), k), k, QR);
	Z=B_vec_remove`/B_SE;
	B_power= 1-probnorm((probit (1-(&alpha/2)))-Z)+ probnorm((-probit (1-(&alpha/2)))-Z);
end;		
if &model= 'random' then do;
 w=1/var_vec;		 
 w_matrix=diag(w);
 w_2=1/(var_vec##2);
 w2_matrix=diag(w_2);
 w_sum=sum(w);
 a=w_sum-(trace((inv(design`*w_matrix*design))*design`*w2_matrix*design));
 QE=es`*(w_matrix-w_matrix*design*(inv(design`*w_matrix*design))*design`*w_matrix)*es;
 tau2=(QE-(n-k-1))/a;
	REV= var_vec + tau2;
	W_new=1/REV;
	W_new_matrix=diag(W_new);
	B_cov= design`*W_new_matrix*design;
	B_cov=inv(B_cov);
	submatrix=B_cov[(k-1):(k+1),(k-1):(k+1)];
	B_SE=sqrt(vecdiag(submatrix));
	B_cov_inv= inv(submatrix);
	QR= B_vec_remove*B_cov_inv*B_vec_remove`;
	QR_power= 1-probchi (cinv ((1-&alpha), k), k, QR);
	Z=B_vec_remove`/B_SE;
	B_power= 1-probnorm((probit (1-(&alpha/2)))-Z)+ probnorm((-probit (1-(&alpha/2)))-Z);
	end;
end;

* +--------------------------------------+
     Test of contrast = 0
   +--------------------------------------+;
if &test= 'contrast' then do;
	if &raw_data= 'no' then do;
		use &Dataset;
 		read all var{&weight} into weight;
 		read all var{&T} into T;
 		read all var{&n1 &n2} into n_vec;
 		read all var{&n1} into n1_vec;
 		read all var{&n2} into n2_vec;
 		read all var{&k} into k_vec;
 		n_groups = nrow(weight);
		common_var_vec = J(n_groups,1,0);
 		group_var = J(n_groups,1,0);
 		group_mn = J(n_groups,1,0); 
 		C2 = J(n_groups,1,0); 

		do groups = 1 to n_groups;
 			C2[groups,1] = weight[groups,1]##2;
			if &eff_type = 'z' then do;
				group_mn[groups,1] = .5 * log((1+T[groups,1]) / (1-T[groups,1]));
				common_var_vec[groups,1] = 1 / (n_vec[groups,1] - 3);
			end;
			if &eff_type=  'r' then do;
				group_mn[groups,1] = T[groups,1];
				common_var_vec[groups,1]=(((1-T[groups,1]**2)**2)/(n_vec[groups,1]-1));
			end; 
			if &eff_type = 'd' then do;
				group_mn[groups,1] = T[groups,1];
				common_var_vec[groups,1] = ((n_vec[groups,1]+n_vec[groups,2])/(n_vec[groups,1]*n_vec[groups,2]))+((T[groups,1]**2)/(2*(n_vec[groups,1]+n_vec[groups,2])));
			end;
			if &eff_type = 'or' then do;
				group_mn[groups,1] = T[groups,1];
				common_var_vec[groups,1] = n_vec[groups,1];
			end;
 		end;
		
		if &model= 'fixed' then do;
			group_var= (common_var_vec/k_vec);
		end;
		if &model= 'random' then do;
				if &tau2^=99 then do;
					total=common_var_vec+&tau2;
					group_var= (total/ k_vec);
				end;
				if &tau2=99 then do;
					total=common_var_vec+(&heterogeneity#common_var_vec);
					group_var= (total/ k_vec);
				end;
		end;
 		G = weight`*group_mn; 
 		Var_G = C2`*group_var; 
 		SE=sqrt(Var_G);
 		Z= G/SE;
 		power= 1-probnorm((probit (1-(&alpha/2)))-Z)+probnorm((-probit (1-(&alpha/2)))-Z);
	end;
	if &raw_data= 'yes' then do;
		use &Dataset;
		read all var{&es} into es; 
		read all var{&v} into var;
		read all var{&weight} into weight; 
 		read all var{&T} into T;
		k= ncol(var);
		n= nrow(var);
		T= remove (T,((k+1):n));
		weight= remove (weight,((k+1):n));
		C2=weight##2;
		inv_var=1/var;
		if &model= 'fixed' then do;
			columnstack_inv=J(1, 2, 0);
			do col=1 to k;
				do row=1 to n;
					if inv_var[row,col] ^=. then do;
						columnstack_inv=columnstack_inv//(col||inv_var[row, col]);
					end;
				end;
			end;
			columnstack_inv=columnstack_inv[2:nrow(columnstack_inv),];
			group_var=J(k,1,0);
				do row=1 to nrow(columnstack_inv);
					group_var[columnstack_inv[row,1],]=group_var[columnstack_inv[row,1],]+columnstack_inv[row,2];
				end;
			group_var=1/group_var;
			G = weight*T`; 
			Var_G = C2*group_var; 
 			SE=sqrt(Var_G);
 			Z= G/SE;
 			power= 1-probnorm((probit (1-(&alpha/2)))-Z)+probnorm((-probit (1-(&alpha/2)))-Z);
			weight=weight`;
			T=T`;
		end;
		If &model= 'random' then do;
			columnstack=J(1, 1, 0);
			do col=1 to k;
				do row=1 to n;
					if es[row,col] ^=. then do;
						columnstack=columnstack//es[row, col];
					end;
				end;
			end;
		columnstack=columnstack[2:nrow(columnstack),];
		columnstack_inv=J(1, 2, 0);
			do col=1 to k;
				do row=1 to n;
					if inv_var[row,col] ^=. then do;
						columnstack_inv=columnstack_inv//(col||inv_var[row, col]);
					end;
				end;
			end;
		columnstack_inv=columnstack_inv[2:nrow(columnstack_inv),];
		columnstack_merge=columnstack_inv||columnstack;
		group_es_w=J(k,1,0);
		group_weight=J(k,1,0);
		group_mean=J(k,1,0);
			do row=1 to nrow(columnstack_merge);
				group_es_w[columnstack_merge[row,1],]=group_es_w[columnstack_merge[row,1],]+(columnstack_merge[row,2]#columnstack_merge[row,3]);
				group_weight[columnstack_merge[row,1],]=group_weight[columnstack_merge[row,1],]+columnstack_merge[row,2];
			end;
		mes=(group_es_w/group_weight)`;
		ones= J(n,1,1);
		m_matrix=ones@mes;
		weight_sq_difference=inv_var#((m_matrix-es)##2);
		Q=weight_sq_difference[+,];
		sum_Q=Q[,+];
		sum_w=inv_var[+,];
		sum_ww=(inv_var##2)[+,];
		a= sum_w-(sum_ww/sum_w);
		sum_a= a[,+];	
		counter=J(n,k,1);
			do col=1 to k;
				do row=1 to n;
					if es[row,col] ^=. then do;
						counter[row,col]=1;
					end;
					if es[row,col] =. then do;
						counter[row,col]=0;
					end;
				end;
			end;
		df= counter[+,]-1;
		sum_df= df[,+];
		tau2=(sum_Q-sum_df)/sum_a;
		if tau2< 0 then do;
			new_tau2=0;
		end;
		if tau2>= 0 then do;
			new_tau2=tau2;
		end;
		random_var= var+new_tau2;
		random_inv_var=1/random_var;
		group_weights=random_inv_var[+,];
		group_var= 1/group_weights;

		G = weight*T`; 
		Var_G = C2*group_var`; 
 		SE=sqrt(Var_G);
 		Z= G/SE;
 		power= 1-probnorm((probit (1-(&alpha/2)))-Z)+probnorm((-probit (1-(&alpha/2)))-Z);
		t_format=t`;
		k_format=(df+1)`;
		weight_format=weight`;
		new_tau2_formatted= new_tau2`; 
		group_var_formatted=group_var`;
		end;
	end;
end;
* +--------------------------------------+
     Test of categorical moderators (QB)=0
   +--------------------------------------+;
if &test= 'QB' then do;
	if &raw_data= 'no' then do;
 		use &Dataset;
 		read all var{&T} into T;
 		read all var{&n1 &n2} into n_vec;
 		read all var{&n1} into n1_vec;
 		read all var{&n2} into n2_vec;
 		read all var{&k} into k_vec;
 		n_groups = nrow(T);
 		common_var_vec = J(n_groups,1,0);
 		group_mn = J(n_groups,1,0); 
 		ones = J(n_groups,1,1);
 		do groups = 1 to n_groups;
			if &eff_type = 'z' then do;
				group_mn[groups,1] = .5 * log((1+T[groups,1]) / (1-T[groups,1]));
				common_var_vec[groups,1] = 1 / (n_vec[groups,1] - 3);	
			end;
			if &eff_type=  'r' then do;
				group_mn[groups,1] = T[groups,1];
				common_var_vec[groups,1]=(((1-T[groups,1]**2)**2)/(n_vec[groups,1]-1));
			end; 
			if &eff_type = 'd' then do;
				group_mn[groups,1] = T[groups,1];
				common_var_vec[groups,1] = ((n_vec[groups,1]+n_vec[groups,2])/(n_vec[groups,1]*n_vec[groups,2]))+((T[groups,1]**2)/(2*(n_vec[groups,1]+n_vec[groups,2])));
			end;
			if &eff_type = 'or' then do;
				group_mn[groups,1] = T[groups,1];
				common_var_vec[groups,1] = n_vec[groups,1];
			end;
 		end;
		if &model= 'fixed' then do;
			total=common_var_vec;
		end;
		if &model= 'random' then do;
			if &tau2^=99 then do;
				total=common_var_vec+&tau2;
			end;
			if &tau2=99 then do;
				total=common_var_vec+(&heterogeneity#common_var_vec);
			end;
		end;
		group_var= (total/ k_vec);
		grand_mn= ((group_mn/group_var)`*ones)/((ones/group_var)`*ones);
		grand_mn_vec = ones*grand_mn;
		noncentrality = (((group_mn-grand_mn_vec)##2)/(group_var))`*ones;
		power= 1-probchi (cinv ((1-&alpha), (n_groups-1)), (n_groups-1), noncentrality);
	end;
	if &raw_data= 'yes' then do;
		use &Dataset;
		read all var{&es} into es; 
		read all var{&v} into var;
 		read all var{&T} into T;
		k= ncol(var);
		n= nrow(var);
		T= remove (T,((k+1):n));
		T_format=T`;
		n_groups=nrow (T_format);
		inv_var=1/var;
		if &model='fixed' then do;
			columnstack_inv=J(1, 2, 0);
				do col=1 to k;
					do row=1 to n;
						if inv_var[row,col] ^=. then do;
							columnstack_inv=columnstack_inv//(col||inv_var[row, col]);
						end;
					end;
				end;
			columnstack_inv=columnstack_inv[2:nrow(columnstack_inv),];
			group_var=J(k,1,0);
					do row=1 to nrow(columnstack_inv);
						group_var[columnstack_inv[row,1],]=group_var[columnstack_inv[row,1],]+columnstack_inv[row,2];
					end;
			group_var=1/group_var;
			inv_group_var=1/group_var;
			sum_inv_group_var=inv_group_var [+,];
			numerator=t_format/group_var;
			numerator_sum= numerator[+,];
			grand_mn= numerator_sum/sum_inv_group_var;
			ones= J(n_groups,1,1);
			grand_mn_vec = ones*grand_mn;
			noncentrality = (((t_format-grand_mn_vec)##2)/(group_var))`*ones;
			power= 1-probchi (cinv ((1-&alpha), (n_groups-1)), (n_groups-1), noncentrality);
			counter=J(n,k,1);
			do col=1 to k;
				do row=1 to n;
					if es[row,col] ^=. then do;
						counter[row,col]=1;
					end;
					if es[row,col] =. then do;
						counter[row,col]=0;
					end;
				end;
			end;
			k_counter=counter[+,];
			k_counter=k_counter`;
		end;

	
		if &model='random' then do;	
			columnstack=J(1, 1, 0);
			do col=1 to k;
				do row=1 to n;
					if es[row,col] ^=. then do;
						columnstack=columnstack//es[row, col];
					end;
				end;
			end;
		columnstack=columnstack[2:nrow(columnstack),];
		columnstack_inv=J(1, 2, 0);
			do col=1 to k;
				do row=1 to n;
					if inv_var[row,col] ^=. then do;
						columnstack_inv=columnstack_inv//(col||inv_var[row, col]);
					end;
				end;
			end;
		columnstack_inv=columnstack_inv[2:nrow(columnstack_inv),];
		columnstack_merge=columnstack_inv||columnstack;
		group_es_w=J(k,1,0);
		group_weight=J(k,1,0);
		group_mean=J(k,1,0);
			do row=1 to nrow(columnstack_merge);
				group_es_w[columnstack_merge[row,1],]=group_es_w[columnstack_merge[row,1],]+(columnstack_merge[row,2]#columnstack_merge[row,3]);
				group_weight[columnstack_merge[row,1],]=group_weight[columnstack_merge[row,1],]+columnstack_merge[row,2];
			end;
		mes=(group_es_w/group_weight)`;
		ones= J(n,1,1);
		m_matrix=ones@mes;
		weight_sq_difference=inv_var#((m_matrix-es)##2);
		Q=weight_sq_difference[+,];
		sum_Q=Q[,+];
		sum_w=inv_var[+,];
		sum_ww=(inv_var##2)[+,];
		a= sum_w-(sum_ww/sum_w);
		sum_a= a[,+];	
		counter=J(n,k,1);
			do col=1 to k;
				do row=1 to n;
					if es[row,col] ^=. then do;
						counter[row,col]=1;
					end;
					if es[row,col] =. then do;
						counter[row,col]=0;
					end;
				end;
			end;
		df= counter[+,]-1;
		sum_df= df[,+];
		tau2=(sum_Q-sum_df)/sum_a;
		if tau2< 0 then do;
			new_tau2=0; 
		end;
		if tau2>= 0 then do;
			new_tau2=tau2;
		end;
		random_var= var+new_tau2;
		random_inv_var=1/random_var;
		group_weights=random_inv_var[+,];
		group_var= 1/group_weights;
		inv_group_var=1/group_var;
		sum_inv_group_var=inv_group_var [,+];
		numerator=t/group_var;
		numerator_sum= numerator[,+]; 
		grand_mn= numerator_sum/sum_inv_group_var;
		ones= J(n_groups,1,1);
		grand_mn_vec = ones*grand_mn;
		k_format=(df+1)`; 
		group_var_formatted=group_var`;
		noncentrality = (((t_format-grand_mn_vec)##2)/(group_var_formatted))`*ones;
		power= 1-probchi (cinv ((1-&alpha), (n_groups-1)), (n_groups-1), noncentrality);
		end;
	end;
end;	
* +----------------------------------------------------------+
   Test of heterogeneity in effect sizes in Regression (QE)= 0
  +----------------------------------------------------------+;
if &test= 'QE' then do;
	if &raw_data='no' then do;
		if &eff_type = 'z' then do;
			z = .5 * log((1+&T) / (1-&T));
			common_var = 1 / (&n1 - 3);
			var_overall = common_var / &k;
		end;
		if &eff_type=  'r' then do;
			common_var=(((1-&T**2)**2)/(&n1-1));
			var_overall = common_var / &k;
		end; 
		if &eff_type = 'd' then do;
			common_var = ((&n1+&n2)/(&n1*&n2))+((&T**2)/(2*(&n1+&n2)));
			var_overall = common_var / &k;
		end;
		if &eff_type = 'or' then do;
			common_var=&n1;
			var_overall = common_var / &k;
		end;
		if &model= 'fixed' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (&k-&p-1)), (&k-&p-1), (&tau2/common_var)*(&k-&p-1));
			end;
			if &tau2 = 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (&k-&p-1)), (&k-&p-1), &heterogeneity*(&k-&p-1));
			end;
		end;
		if &model= 'random' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (((cinv ((1-&alpha), (&k-&p))*common_var)/(common_var+ &tau2)), (&k-&p-1), 0);
			end;
			if &tau2 =99 then do;
				power= 1-probchi (((cinv ((1-&alpha), (&k-&p))*common_var)/(common_var+ (common_var*&heterogeneity))), (&k-&p-1), 0);	
			end;
		end;
	end;
	if &raw_data='yes' then do;
		use &Dataset;
		read all var{&x} into x_x;
		read all var{&v}  into v; 
		read all var{&B}  into B;
 		read all var{&es}  into es;
		p= ncol(x_x);
 		k= nrow(x_x);
		w=1/v;
		var_overall=1/sum(w);
		average_var=var_overall*k;
		if &model= 'fixed' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (k-p-1)), k-p-1, (&tau2/average_var)*(k-p-1));
			end;
			if &tau2 = 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (k-p-1)), k-p-1, ((&heterogeneity*average_var)/average_var)*(k-p-1));
			end;
		end;
		if &model= 'random' then do;
		x = j(k,1,1) || x_x ;
		B_vec_remove= remove (B,((p+1):k));			 
		w_matrix=diag(w);
		w_2=1/(v##2);
		w2_matrix=diag(w_2);
		w_sum=sum(w);
		a=w_sum-(trace((inv(x`*w_matrix*x))*x`*w2_matrix*x));
		M=w_matrix*x*(inv(x`*w_matrix*x))*x`*w_matrix;	
		v_random=v+&tau2;
		w_random=1/v_random;
		w_random_matrix=diag(w_random);
		w2_random=1/(v_random##2);
		w2_random_matrix=diag(w2_random);
		mean=(a*&tau2)+(k-p-1);
		/*variance=2*((trace(w2_matrix*w2_random_matrix))-(2*(trace(M*w_matrix*w2_random_matrix)))+(trace(M*w_random_matrix*M*w_random_matrix)));
*/	
		variance= 128.43;
		r=variance/(2*mean);
		s=(2*(mean##2))/variance;
		power=1-probchi (cinv ((1-&alpha), (k-p-1))/r, s, 0);
		end;
	end;
end;
* +---------------------------------------------------------------+
   Test of heterogeneity in effect sizes across all levels (QW)= 0
  +---------------------------------------------------------------+;
if &test= 'QW' then do;
	if &raw_data= 'no' then do;
		if &eff_type = 'z' then do;
			z = .5 * log((1+&T) / (1-&T));
			common_var = 1 / (&n1 - 3);
			var_overall = common_var / &k;
		end;
		if &eff_type=  'r' then do;
			common_var=(((1-&T**2)**2)/(&n1-1));
			var_overall = common_var / &k;
		end; 
		if &eff_type = 'd' then do;
			common_var = ((&n1+&n2)/(&n1*&n2))+((&T**2)/(2*(&n1+&n2)));
			var_overall = common_var / &k;
		end;
		if &eff_type = 'or' then do;
			common_var=&n1;
			var_overall = common_var / &k;
		end;
		if &model= 'fixed' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (&k-&p)), (&k-&p), (&tau2/common_var)*(&k-&p));
			end;
			if &tau2 = 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (&k-&p)), (&k-&p), &heterogeneity*(&k-&p));
				/*common_var=99;
				var_overall=99;*/
			end;
		end;		
		if &model= 'random' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (((cinv ((1-&alpha), (&k-&p))*common_var)/(common_var+ &tau2)), (&k-&p), 0);
			end;
			if &tau2 =99 then do;
				power= 1-probchi (((cinv ((1-&alpha), (&k-&p))*common_var)/(common_var+ (common_var*&heterogeneity))), (&k-&p), 0);
			end;
		end;
	end;
	if &raw_data= 'yes' then do;
			use &Dataset;
			read all var{&es} into es;
			read all var{&v}  into v; 		 
			inv_var=1/v;
			k_= ncol(es);
			n= nrow(es);
			counter=J(n,k_,1);
				do col=1 to k_;
					do row=1 to n;
						if es[row,col] ^=. then do;
							counter[row,col]=1;
						end;
						if es[row,col] =. then do;
							counter[row,col]=0;
						end;
					end;
				end;
			df= counter[+,]-1;
			sum_df= df[,+];
			k=sum_df+k_;
			sum_w=inv_var[+,];
			sum_w_tot=sum_w[,+];
			var_overall=1/sum_w_tot;
			average_var=var_overall*k;
			sum_ww=(inv_var##2)[+,];
			sum_www=(inv_var##3)[+,];
			a= sum_w-(sum_ww/sum_w);
			sum_a= a[,+];	
		if &model= 'fixed' then do;
			if &tau2 ^= 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (sum_df)), (sum_df), (&tau2/average_var)*(sum_df));
			end;
			if &tau2= 99 then do;
				power= 1-probchi (cinv ((1-&alpha), (sum_df)), (sum_df),((&heterogeneity*average_var)/average_var)*(sum_df));
			end;
		end;
		if &model= 'random' then do;
			if &tau2 ^= 99 then do;
				mean=(sum_a*&tau2)+(sum_df); 
				b=(df)+(2*a*&tau2)+(sum_ww-(2*(sum_www/sum_w))+((sum_ww##2)/(sum_w##2)))*(&tau2##2);
			end;
			if &tau2= 99 then do;
				mean=(sum_a*(&heterogeneity*average_var))+(sum_df); 
				b=(df)+(2*a*(&heterogeneity*average_var))+(sum_ww-(2*(sum_www/sum_w))+((sum_ww##2)/(sum_w##2)))*((&heterogeneity*average_var)##2);
			end;
				b_sum= b[,+];
				variance=2*b_sum;
				r=variance/(2*mean);
				s=(2*(mean##2))/variance;
				power=1-probchi (cinv ((1-&alpha), sum_df)/r, s, 0);
		end;
	end;	
end;

* +--------------------+
    Print Statements 
  +--------------------+;
%if &test= 'Reg' %then %do;
print "--------- Meta-Analysis Power Macro ------" ;
print &model [label= ' '] [rowname = 'Model= '] ;
print n [label= ' '] [rowname = 'Number of Studies = '];
%if &model='fixed'%then %do;
print &tau2[label= ' '] [rowname = 'Random Effect Size Variance=']; 
%end;
%if &model='random'%then %do;
print tau2[label= ' '] [rowname = 'Random Effect Size Variance=']; 
%end;
print &alpha [label= ' '] [rowname = 'Alpha = '];
print QR_power [label= ' '] [rowname = 'Power for Overall Model (QR) = '][format=12.8];
print 'Power for Regression Coefficients (B)';
names = {&x} ;
print B_power [rowname=names] [colname=" "] [label=""][format=12.8];
print "--------------------------------------------" ;
%end;

%if &test= 'QB' %then %do;
%if &raw_data='no'|&raw_data='yes'%then %do;
print "-------------------------------Meta-Analysis Power Macro-------------------------------" ;
print "Test of Categorical Moderators" ;
names = &model//&eff_type;
print names [label= ' '] 
[rowname = {
'Model = ' 
'Effect Size Metric = '}];
print &alpha [label= ' '][rowname = 'Alpha = '];
print n_groups [label= ' '][rowname = 'Number of Categories = '];
%if &raw_data='no'%then %do;
Print '********Raw data not provided********';
Print T [label='ES'] 
 n1_vec [label='n(g1)'] 
 n2_vec [label='n(g2)']
 k_vec [label='K'] 
 common_var_vec [label='Common Variance'] 
 &tau2 [label='Tau-Squared'] 
 &heterogeneity [label='Heterogeneity'] 
 group_var [label='Group Variance']; 
print power [label= ' '][rowname = 'Estimated Power of Test = '];
print "----------------------------------------------------------------------------------------" ;
%end;
%if &raw_data='yes'%then %do;
Print '********Raw data provided********';
%if &model='fixed'%then %do;
 Print T_format [label='ES'] 
 k_counter [label='K']
 group_var [label='Group Variance']; 
print power [label= ' '][rowname = 'Estimated Power of Test = '];
print "----------------------------------------------------------------------------------------" ;
%end;
%if &model='random'%then %do;
Print t_format [label='ES']
k_format [label='K'] 
new_tau2 [label='Tau-Squared']
group_var_formatted [label='Group Variance']; 
print power [label= ' '][rowname = 'Estimated Power of Test = '];
print "--------------------------------------------------------------------------------------" ;
%end;
%end;
%end;
%end;


%if &test= 'contrast' %then %do;
%if &raw_data='no'|&raw_data='yes'%then %do;
print "-------------------------------Meta-Analysis Power Macro-------------------------------" ;
print "Test of Contrasts" ;
names = &model//&eff_type;
print names [label= ' '] 
[rowname = {
'Model = ' 
'Effect Size Metric = '}];
print &alpha [label= ' '][rowname = 'Alpha = '];
%if &raw_data='no'%then %do;
 Print '********Raw data not provided********';
 Print T [label='ES'] 
 weight  [label='weight']
 n1_vec [label='n(g1)'] 
 n2_vec [label='n(g2)']
 k_vec [label='K'] 
 common_var_vec [label='Common Variance'] 
 &tau2 [label='Tau-Squared'] 
 &heterogeneity [label='Heterogeneity'] 
 group_var [label='Group Variance']; 
print power [label= ' '][rowname = 'Estimated Power of Test = '];
print "----------------------------------------------------------------------------------------" ;
%end;
%if &raw_data='yes'%then %do;
print '********Raw data provided********';
%if &model='fixed'%then %do;
Print T [label='ES'] 
weight  [label='weight']
group_var [label='Group Variance']; 
print power [label= ' '][rowname = 'Estimated Power of Test = '];
print "----------------------------------------------------------------------------------------" ;
%end;
%if &model='random'%then %do;
Print t_format [label='ES']
weight_format  [label='weight']
k_format [label='K'] 
new_tau2_formatted [label='Tau-Squared']
group_var_formatted [label='Group Variance']; 
print power [label= ' '][rowname = 'Estimated Power of Test = '];
print "----------------------------------------------------------------------------------------" ;
%end;
%end;
%end;
%end;

%if &test= 'M'%then %do;
print "---------------------Meta-Analysis Power Macro---------------------" ;
names = &test//&model//&eff_type;
print names [label= ' '] 
[rowname = {'Significance Test Type = ' 
'Model = ' 
'Effect Size Metric = '}];
%if &raw_data='no'%then %do;
print '********Raw data not provided********';
values = &T//&n1//&n2//&k//common_var//var_overall//&tau2//&heterogeneity//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Individual Study Sample Size (group 1) = '
'Individual Study Sample Size (group 2) = '
'Number of Studies = '
'Estimated Common Within Study Variance = '
'Estimated Overall Within Study Variance = '
'Random Effects Variance ='
'Heterogeneity Ratio = '
'Alpha = '
'Estimated Power of Test = '}];
print "-------------------------------------------------------------------" ;
%end;
%if &raw_data='yes'%then %do;
%if &model='fixed'%then %do;
print '********Raw data provided********';
values = &T//k//var_overall//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Number of Studies = '
'Overall Within Study Variance = '
'Alpha = '
'Estimated Power of Test = '}];
print"-------------------------------------------------------------------" ;
%end;
%if &model='random'%then %do;
print '********Raw data provided********';
values = &T//k//var_random_overall//tau2//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Number of Studies = '
'Overall Within Study Variance = '
'Random Effects Variance Component='
'Alpha = '
'Estimated Power of Test = '}];
print"-------------------------------------------------------------------" ;
%end;
%end;
%end;

%if &test= 'QE'%then %do;
print "---------------------Meta-Analysis Power Macro---------------------" ;
names = &test//&model//&eff_type;
print names [label= ' '] 
[rowname = {'Significance Test Type = ' 
'Model = ' 
'Effect Size Metric = '}];
%if &raw_data='no'%then %do;
print '********Raw data not provided********';
values = &T//&n1//&n2//&k//common_var//var_overall//&tau2//&heterogeneity//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Individual Study Sample Size (group 1) = '
'Individual Study Sample Size (group 2) = '
'Number of Studies = '
'Estimated Common Within Study Variance = '
'Estimated Overall Within Study Variance = '
'Random Effects Variance ='
'Heterogeneity Ratio = '
'Alpha = '
'Estimated Power of Test = '}];
print "-------------------------------------------------------------------" ;
%end;
%if &raw_data='yes'%then %do;
%if &model='fixed'%then %do;
print '********Raw data provided********';
values = &T//k//average_var//var_overall//&tau2//&heterogeneity//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Number of Studies = '
'Common Variance='
'Overall Within Study Variance = '
'Ranom Effects Variance'
'Heterogeneity'
'Alpha = '
'Estimated Power of Test = '}];
print"-------------------------------------------------------------------" ;
%end;
%if &model='random'%then %do;
print variance;
print '********Raw data provided********';
values = &T//k//&tau2//&heterogeneity//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Number of Studies = '
'Random Effects Variance'
'Heterogeneity Ratio'
'Alpha = '
'Estimated Power of Test = '}];
print"-------------------------------------------------------------------" ;
%end;
%end;
%end;

%if &test= 'QW'%then %do;
print "---------------------Meta-Analysis Power Macro---------------------" ;
names = &test//&model//&eff_type;
print names [label= ' '] 
[rowname = {'Significance Test Type = ' 
'Model = ' 
'Effect Size Metric = '}];
%if &raw_data='no'%then %do;
print '********Raw data not provided********';
values = &T//&n1//&n2//&k//common_var//var_overall//&tau2//&heterogeneity//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Individual Study Sample Size (group 1) = '
'Individual Study Sample Size (group 2) = '
'Number of Studies = '
'Estimated Common Within Study Variance = '
'Estimated Overall Within Study Variance = '
'Random Effects Variance ='
'Heterogeneity Ratio = '
'Alpha = '
'Estimated Power of Test = '}];
print "-------------------------------------------------------------------" ;
%end;
%if &raw_data='yes'%then %do;
%if &model='fixed'%then %do;
print '********Raw data provided********';
values = &T//k//var_overall//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Number of Studies = '
'Overall Within Study Variance = '
'Alpha = '
'Estimated Power of Test = '}];
print"-------------------------------------------------------------------" ;
%end;
%if &model='random'%then %do;
print '********Raw data provided********';
values = &T//k//&tau2//&heterogeneity//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Number of Studies = '
'Random Effects Variance'
'Heterogeneity Ratio'
'Alpha = '
'Estimated Power of Test = '}];
print"-------------------------------------------------------------------" ;
%end;
%end;
%end;

%if &test= 'QT'%then %do;
print "---------------------Meta-Analysis Power Macro---------------------" ;
names = &test//&model//&eff_type;
print names [label= ' '] 
[rowname = {'Significance Test Type = ' 
'Model = ' 
'Effect Size Metric = '}];
%if &raw_data='no'%then %do;
print '********Raw data not provided********';
values = &T//&n1//&n2//&k//common_var//var_overall//&tau2//&heterogeneity//&alpha//power;
print values [label= ' '] [rowname = {
'Population Effect Size = '
'Individual Study Sample Size (group 1) = '
'Individual Study Sample Size (group 2) = '
'Number of Studies = '
'Estimated Common Within Study Variance = '
'Estimated Overall Within Study Variance = '
'Random Effects Variance='
'Heterogeneity Ratio = '
'Alpha = '
'Estimated Power of Test = '}];
print "-------------------------------------------------------------------" ;
%end;
%if &raw_data='yes'%then %do;
%if &model='fixed'%then %do;
print '********Raw data provided********';
values = &tau2//&heterogeneity//k//&alpha//power;
print values [label= ' '] [rowname = {
'Random Effects Variance='
'Heterogeneity Ratio =' 
'Number of Studies = '
'Alpha = '
'Estimated Power of Test = '}];
print"-------------------------------------------------------------------" ;
%end;
%if &model='random'%then %do;
print '********Raw data provided********';
values = &tau2//&heterogeneity//k//&alpha//power;
print values [label= ' '] [rowname = {
'Random Effects Variance='
'Heterogeneity Ratio =' 
'Number of Studies = '
'Alpha = '
'Estimated Power of Test = '}];
print"-------------------------------------------------------------------" ;
%end;
%end;
%end;
quit;
%mend metapower;

******************************************M*****************************************************;

%metapower (test= 'M', model= 'fixed', raw_data= 'no' , alpha= .05, tau2= 99, heterogeneity = 99, n1= 12, n2= 36, k= 18,
eff_type= 'd' , T= .20, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*H&P (2001), p.207;
data hedges1;
input v es;
cards;
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
.111 .20
;

%metapower(test= 'M', model= 'fixed', raw_data= 'yes', alpha= .05, tau2= 99, heterogeneity= 99, n1= 99, n2= 99, k= 99,
eff_type= 'd' , T= .20, Dataset= hedges1, B= NA, v= v, x= NA, es= es, p= NA, weight=NA );run;*same as previous statement except here data are available;


%metapower(test= 'M', model= 'fixed', raw_data= 'no', alpha= .05, tau2= 99, heterogeneity= 99, n1= 25, n2= 0, k= 10,
eff_type= 'z' , T= .10, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*H&P (2001), p.208;

data hedges2;
input v es;
cards;
.045 .10
.045 .10
.045 .10
.045 .10
.045 .10
.045 .10
.045 .10
.045 .10
.045 .10
.045 .10
;

%metapower(test= 'M', model= 'fixed', raw_data= 'yes', alpha= .05, tau2= 99, heterogeneity= 99, n1= 99, n2=99, k= 99,
eff_type= 'z' , T= .10, Dataset= hedges2, B= NA, v= v, x= NA, es= es, p= NA, weight=NA );run;*same as previous statement except here data are available;

%metapower(test= 'M', model= 'random', raw_data= 'no', alpha= .05, tau2= .037, heterogeneity= 99, n1= 12, n2=36, k= 18,
eff_type= 'd' , T= .20, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*H&P (2001), p.213;

%metapower(test= 'M', model= 'random', raw_data= 'no', alpha= .05, tau2= 99, heterogeneity= .33, n1= 12, n2=36, k= 18,
eff_type= 'd' , T= .20, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*Same as previous statement except heterogeneity ratio is used instead of tau-squared;


data random_data1;
input es v;
cards;
0.42364893	0.005
0.447692024	0.005714286
0.484700279	0.004
0.693147181	0.004
0.618381314	0.005
0.775298706	0.004444444
;

%metapower(test= 'M', model= 'fixed', raw_data= 'yes', alpha= .05, tau2= 99, heterogeneity= 99, n1= 99, n2=99, k= 99,
eff_type= 'z' , T= .05, Dataset= random_data1, B= NA, v= v, x= NA, es= es, p= NA, weight=NA );run;*Brannick website example;
%metapower(test= 'M', model= 'random', raw_data= 'yes', alpha= .05, tau2= 99, heterogeneity= 99, n1= 99, n2=99, k= 99,
eff_type= 'z' , T= .05, Dataset= random_data1, B= NA, v= v, x= NA, es= es, p= NA, weight=NA );run;*Brannick website example;

******************************************QT*********************************************************;

%metapower(test= 'QT', model= 'fixed', raw_data= 'no', alpha= .05, tau2= 99, heterogeneity= .67, n1= 25, n2=0, k= 10,
eff_type= 'z' , T= .10, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*H&P (2001), p.209-210;
%metapower(test= 'QT', model= 'fixed', raw_data= 'no', alpha= .05, tau2= .03045, heterogeneity= 99, n1= 25, n2=0, k= 10,
eff_type= 'z' , T= .10, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*Same as previous statement except tau-squared is used instead of heterogeneity ratio;
%metapower(test= 'QT', model= 'fixed', raw_data= 'yes', alpha= .05, tau2= 99, heterogeneity= .67, n1= 99, n2=99, k= 99,
eff_type= 'z' , T= 99, Dataset= random_data1, B= NA, v= v, x= NA, es= es, p= NA, weight=NA );run;*Brannick website example;
%metapower(test= 'QT', model= 'random', raw_data= 'no', alpha= .05, tau2= .015, heterogeneity= 99, n1= 25, n2=0, k= 10,
eff_type= 'z' , T= .10, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*H&P (2001), p.214;
%metapower(test= 'QT', model= 'random', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=.33, n1= 25, n2=0, k= 10,
eff_type= 'z' , T= .10, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*Same as previous statement except heterogeneity ratio is used instead of tau-squared;

%metapower(test= 'QT', model= 'random', raw_data= 'yes', alpha= .05, tau2=.005, heterogeneity=99, n1= 99, n2=99, k= 99,
eff_type= 'z' , T= .10, Dataset=random_data1, B= NA, v= v, x= NA, es= es, p= NA, weight=NA );run;*Brannick website data with random-effects w/ data set;

%metapower(test= 'QT', model= 'random', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=1.0684, n1= 99, n2=99, k= 99,
eff_type= 'z' , T= .10, Dataset=random_data1, B= NA, v= v, x= NA, es= es, p= NA, weight=NA );run;*Same as previous statement except heterogeneity ratio is used instead of tau-squared;


******************************************Contrasts*********************************************************;
data hedges3;
input weight T nn1 nn2 kk;
cards;
1  .49 30  30  8
-1 .29 107 107 29  
;


%metapower(test= 'contrast', model= 'fixed', raw_data= 'no', alpha= .10, tau2=99, heterogeneity=99, n1= nn1, n2=nn2, k= kk,
eff_type= 'd' , T= T, Dataset=hedges3, B= NA, v= NA, x= NA, es= NA, p= NA, weight=weight );run;

%metapower(test= 'contrast', model= 'random', raw_data= 'no', alpha= .10, tau2=.10, heterogeneity=99, n1= nn1, n2=nn2, k= kk,
eff_type= 'd' , T= T, Dataset=hedges3, B= NA, v= NA, x= NA, es= NA, p= NA, weight=weight );run;

data hedges4;
input es1 v1 es2 v2 weight T ;
cards;
.49 .069 .29 .02 1 .49 
.49 .069 .29 .02 -1 .29   
.49 .069 .29 .02 . .
.49 .069 .29 .02 . .
.49 .069 .29 .02 . .
.49 .069 .29 .02 . .
.49 .069 .29 .02 . .
.49 .069 .29 .02 . .
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
.   .    .29 .02 . . 
;


%metapower(test= 'contrast', model= 'fixed', raw_data= 'yes', alpha= .10, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd' , T= T, Dataset=hedges4, B= NA, v= v1 v2, x= NA, es= es1 es2, p= NA, weight=weight );run;*same as previous statement except here data are available;



Data hedges5;
input es1 v1 es2 v2 es3 v3 T weight ;
cards;
.72  .03  .27  .06 .51 .03 .25  1
.63  .04  -.11 .05 .25 .01 .00  -.5
.73  .04  .00  .17 .38 .01 .00  -.5
.07  .11  .53  .07 .60 .02 .    .
1.19 .15  .14  .17 .91 .06 .    .  
.47  .14  -.07 .17 .36 .06 .    .
.    .    .39  .03 .12 .11 .    .
.    .    .16  .11 .03 .07 .    .
.    .    .53  .04 .20 .04 .    .
.    .    2.27 .08 .60 .11 .    . 
.    .    .    .   .21 .13 .    . 
.    .    .    .   .24 .10 .    .
.    .    .    .   .50 .10 .    . 
.    .    .    .   .76 .09 .    .
;
%metapower(test= 'contrast', model= 'fixed', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= T, Dataset=hedges5, B= NA, v= v1 v2 v3, x= NA, es= es1 es2 es3, p= NA, weight=weight );run;* H&P (2004), p. 432 (Note difference of .64 and .66 is due to greater computational accuracy);

%metapower(test= 'contrast', model= 'random', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= T, Dataset=hedges5, B= NA, v= v1 v2 v3, x= NA, es= es1 es2 es3, p= NA, weight=weight );run;* H&P (2004), p. 437 ;



******************************************QB*********************************************************;

data hedges6;
input T nn1 nn2 kk;
cards;
.70  38 38 6
.45  30 30 10
;

%metapower(test= 'QB', model= 'fixed', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=99, n1= nn1, n2=nn2, k=kk,
eff_type= 'd', T= T, Dataset=hedges6, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;* if variances unknown from H&P (2004, p.430) but one had good estimates of v based on average study ns and k, one might set up power calculations like this;



Data hedges7;
input es1 v1 es2 v2 T ;
cards;
.72  .03  .27  .06 .00
.63  .04  -.11 .05 .25
.73  .04  .00  .17 .
.07  .11  .53  .07 .
1.19 .15  .14  .17 .  
.47  .14  -.07 .17 .
.    .    .39  .03 .
.    .    .16  .11 .
.    .    .53  .04 .
.    .    2.27 .08 .
;
%metapower(test= 'QB', model= 'fixed', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= T, Dataset=hedges7, B= NA, v= v1 v2, x= NA, es= es1 es2, p= NA, weight=NA );run;* H&P (2004), p. 430 ;

Data hedges8;
input es1 v1 es2 v2 es3 v3 T;
cards;
.72  .03  .27  .06 .51 .03 .25 
.63  .04  -.11 .05 .25 .01 .00 
.73  .04  .00  .17 .38 .01 .125
.07  .11  .53  .07 .60 .02 .
1.19 .15  .14  .17 .91 .06 .  
.47  .14  -.07 .17 .36 .06 .
.    .    .39  .03 .12 .11 .
.    .    .16  .11 .03 .07 .
.    .    .53  .04 .20 .04 .
.    .    2.27 .08 .60 .11 .
.    .    .    .   .21 .13 .
.    .    .    .   .24 .10 .
.    .    .    .   .50 .10 .
.    .    .    .   .76 .09 .
;

%metapower(test= 'QB', model= 'fixed', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= T, Dataset=hedges8, B= NA, v= v1 v2 v3, x= NA, es= es1 es2 es3, p= NA, weight=NA );run;* H&P (2004), p. 430-431;

Data hedges9;
input T nn1 nn2 kk;
cards;
.25   25 25 10
.125  25 25 10
.00   25 25 10
;

%metapower(test= 'QB', model= 'random', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=.33, n1= nn1, n2=nn2, k=kk,
eff_type= 'd', T= T, Dataset=hedges9, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;* H&P (2004), p. 436;

Data hedges10;
input es1 v1 es2 v2 es3 v3 T;
cards;
.72  .03  .27  .06 .51 .03 .25 
.63  .04  -.11 .05 .25 .01 .00 
.73  .04  .00  .17 .38 .01 .125
.07  .11  .53  .07 .60 .02 .
1.19 .15  .14  .17 .91 .06 .  
.47  .14  -.07 .17 .36 .06 .
.    .    .39  .03 .12 .11 .
.    .    .16  .11 .03 .07 .
.    .    .53  .04 .20 .04 .
.    .    2.27 .08 .60 .11 .
.    .    .    .   .21 .13 .
.    .    .    .   .24 .10 .
.    .    .    .   .50 .10 .
.    .    .    .   .76 .09 .
;

%metapower(test= 'QB', model= 'random', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= T, Dataset=hedges10, B= NA, v= v1 v2 v3, x= NA, es= es1 es2 es3, p= NA, weight=NA );run;* hypothetical example;

Data hedges11;
input es1 v1 es2 v2 es3 v3 T weight ;
cards;
.72  .03  .27  .06 .51 .03 .25  1
.63  .04  -.11 .05 .25 .01 .00  -1
.73  .04  .00  .17 .38 .01 .00  0
.07  .11  .53  .07 .60 .02 .    .
1.19 .15  .14  .17 .91 .06 .    .  
.47  .14  -.07 .17 .36 .06 .    .
.    .    .39  .03 .12 .11 .    .
.    .    .16  .11 .03 .07 .    .
.    .    .53  .04 .20 .04 .    .
.    .    2.27 .08 .60 .11 .    . 
.    .    .    .   .21 .13 .    . 
.    .    .    .   .24 .10 .    .
.    .    .    .   .50 .10 .    . 
.    .    .    .   .76 .09 .    .
;
%metapower(test= 'contrast', model= 'random', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= T, Dataset=hedges11, B= NA, v= v1 v2 v3, x= NA, es= es1 es2 es3, p= NA, weight=weight );run;* H&P (2004), p. 435-436 ;

******************************************QW*********************************************************;


%metapower(test= 'QW', model= 'fixed', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=.33, n1= 99, n2=99, k=30,
eff_type= 'd', T= 99, Dataset=NA, B= NA, v= NA, x= NA, es=NA, p= 3, weight=NA );run;* H&P (2004), p. 433 ;
%metapower(test= 'QW', model= 'fixed', raw_data= 'no', alpha= .05, tau2=.034, heterogeneity=99, n1= 20, n2=20, k=30,
eff_type= 'd', T= .50, Dataset=NA, B= NA, v= NA, x= NA, es=NA, p= 3, weight=NA );run;*similar to above statement except tau2 and average within study variance specified; 
%metapower(test= 'QW', model= 'random', raw_data= 'no', alpha= .05, tau2=.027, heterogeneity=99, n1= 25, n2=25, k=30,
eff_type= 'd', T= .50, Dataset=NA, B= NA, v= NA, x= NA, es=NA, p= 3, weight=NA );run;*H&P (2004), p. 439;
%metapower(test= 'QW', model= 'random', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=.33, n1= 25, n2=25, k=30,
eff_type= 'd', T= .50, Dataset=NA, B= NA, v= NA, x= NA, es=NA, p= 3, weight=NA );run;*H&P (2004), p. 439;

data hedges12;
input es1 v1 es2 v2 es3 v3;
cards;
.72  .03  .27  .06 .51 .03 
.63  .04  -.11 .05 .25 .01 
.73  .04  .00  .17 .38 .01 
.07  .11  .53  .07 .60 .02 
1.19 .15  .14  .17 .91 .06 
.47  .14  -.07 .17 .36 .06 
.    .    .39  .03 .12 .11 
.    .    .16  .11 .03 .07 
.    .    .53  .04 .20 .04 
.    .    2.27 .08 .60 .11 
.    .    .    .   .21 .13 
.    .    .    .   .24 .10 
.    .    .    .   .50 .10  
.    .    .    .   .76 .09 
;

%metapower(test= 'QW', model= 'random', raw_data= 'yes', alpha= .05, tau2=.027, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= 99, Dataset=Hedges12, B= NA, v= v1 v2 v3, x= NA, es= es1 es2 es3, p= 3, weight=NA );run;*H&P (2004), p. 438-439;

%metapower(test= 'QW', model= 'random', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=.60, n1= 99, n2=99, k=99,
eff_type= 'd', T= 99, Dataset=Hedges12, B= NA, v= v1 v2 v3, x= NA, es= es1 es2 es3, p= 3, weight=NA );run;*same as above except heterogeneity used; 
******************************************QE*********************************************************;
%metapower(test= 'QE', model= 'fixed', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=.33, n1= 25, n2=25, k=19,
eff_type= 'r', T= .10, Dataset=NA, B= NA, v= NA, x= NA, es= NA, p=3, weight=NA );run;*H&P (2004)pg.441-442;

%metapower(test= 'QE', model= 'fixed', raw_data= 'no', alpha= .05, tau2=.0136, heterogeneity=99, n1= 25, n2=25, k=19,
eff_type= 'r', T= .10, Dataset=NA, B= NA, v= NA, x= NA, es= NA, p=3, weight=NA );run;*same as above except tau-squared specified;

%metapower(test= 'QE', model= 'random', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=.33, n1= 25, n2=25, k=19,
eff_type= 'r', T= .10, Dataset=NA, B= NA, v= NA, x= NA, es= NA, p=3, weight=NA );run;

%metapower(test= 'QE', model= 'random', raw_data= 'no', alpha= .05, tau2=.0136, heterogeneity=99, n1= 25, n2=25, k=19,
eff_type= 'r', T= .10, Dataset=NA, B= NA, v= NA, x= NA, es= NA, p=3, weight=NA );run;

data hedges13;
input B es v x1 x2 x3;
cards;

.250 -.17  .03 0 1 0 
.250 .38   .03 0 0 0 
0.00 .23   .01 1 0 0 
.    .02   .09 1 1 0 
.    1.92  .07 1 0 0 
.	 1.63  .08 1 1 0 
.	 .56   .06 1 1 1 
.	 -1.11 .19 1 0 1
.	 .20   .17 1 0 1 
.	 -.08  .17 1 0 1 
.	 .53   .16 0 1 0 
.	 .93   .09 1 1 1 
.	 .56   .04 0 1 0 
.	 .08   .06 0 1 0 
.	 .52   .06 0 1 1 
.	 .07   .11 1 1 0
.    2.94  .13 1 1 1 
.	 1.30  .16 1 1 0 
.	 .04   .14 0 1 1 
;

%metapower(test= 'QE', model= 'fixed', raw_data= 'yes', alpha= .05, tau2=.067, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'r', T= .10, Dataset=hedges13, B= B, v= v, x= x1 x2 x3, es= es, p=99, weight=NA );run;*hypothetical example w/ data fixed;
%metapower(test= 'QE', model= 'fixed', raw_data= 'yes', alpha= .05, tau2= 99, heterogeneity=1.21, n1= 99, n2=99, k=99,
eff_type= 'r', T= .10, Dataset=hedges13, B= B, v= v, x= x1 x2 x3, es= es, p=99, weight=NA );run;*hypothetical example w/ data fixed;
%metapower(test= 'QE', model= 'random', raw_data= 'yes', alpha= .05, tau2= .067, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'r', T= .10, Dataset=hedges13, B= B, v= v, x= x1 x2 x3, es= es, p=99, weight=NA );run;*H&P (2004) p.444;

******************************************Fixed Parameters in Regression*********************************************************;

%metapower(test= 'Reg', model= 'fixed', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= 99, Dataset=Hedges13, B= B, v= v, x= x1 x2 x3, es= es, p= NA, weight=NA );run;*H&P (2004) 440-441;

%metapower(test= 'Reg', model= 'random', raw_data= 'yes', alpha= .05, tau2=99, heterogeneity=99, n1= 99, n2=99, k=99,
eff_type= 'd', T= 99, Dataset=Hedges13, B= B, v= v, x= x1 x2 x3, es= es, p= NA, weight=NA );run;


********************************************Odds Ratio Examples*******************************************************************;
%metapower (test= 'M', model= 'random', raw_data= 'no' , alpha= .05, tau2= 99, heterogeneity = 99, n1= .10, n2= 99, k= 18,
eff_type= 'or' , T= -1.045, Dataset= odds, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA);run;
%metapower(test= 'QT', model= 'fixed', raw_data= 'no', alpha= .05, tau2= 99, heterogeneity= .67, n1= .10, n2=0, k= 10,
eff_type= 'or' , T= -1.045, Dataset= NA, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;*H&P (2001), p.209-210;
data odds1;
input weight T nn1 nn2 kk;
cards;
1  -1.045 .10  0  8
-1 -1.045 .10  0  29  
;


%metapower(test= 'contrast', model= 'fixed', raw_data= 'no', alpha= .10, tau2=99, heterogeneity=99, n1= nn1, n2=nn2, k= kk,
eff_type= 'or' , T= T, Dataset=odds1, B= NA, v= NA, x= NA, es= NA, p= NA, weight=weight );run;
data odds2;
input T nn1 nn2 kk;
cards;
-1.045 .10 0 6
-1.045 .10 0 10
;

%metapower(test= 'QB', model= 'fixed', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=99, n1= nn1, n2=nn2, k=kk,
eff_type= 'or', T= T, Dataset=odds2, B= NA, v= NA, x= NA, es= NA, p= NA, weight=NA );run;* if variances unknown from H&P (2004, p.430) but one had good estimates of v based on average study ns and k, one might set up power calculations like this;

%metapower(test= 'QW', model= 'fixed', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=.33, n1= .10, n2=99, k=30,
eff_type= 'or', T= -1.045, Dataset=NA, B= NA, v= NA, x= NA, es=NA, p= 3, weight=NA );run;* H&P (2004), p. 433 ;

%metapower(test= 'QE', model= 'fixed', raw_data= 'no', alpha= .05, tau2=99, heterogeneity=.33, n1= .10, n2=25, k=19,
eff_type= 'or', T= -1.045, Dataset=NA, B= NA, v= NA, x= NA, es= NA, p=3, weight=NA );run;*H&P (2004)pg.441-442;






