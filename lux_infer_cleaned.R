# Lux inference code
# Matlab code by Mudassar Iqbal
# R translation by Dov Stekel
# Funded by BBSRC BB/I001875/1


## model ##########

require(deSolve)

# ===========Constant Parameters===========================================
luxparams = c(
# ============V_max used in velocity equations.================
Vm_Fre = 61.1, # VM*Fre                                                           
Vm_LuxAB = 64.97, #?
Vm_LuxEC = 198.93, #?
# ============Fixed Concentrations ==============================
Conc_F = 88,
Conc_R = 231, 
Conc_NADPH = 560, 
Conc_O2 = 214,
Conc_ATP = 1310,
# ============Fixed Concentrations ==============================
Conc_F = 88, 
Conc_R = 231, 
Conc_NADPH = 560, 
Conc_O2 = 214,
Conc_ATP = 1310,
# For Fre
Fre_kr1 = 1, # Assumed in equilibrium
Fre_k31 =  1.8, #inferred
Fre_k32 = 45.2, # inferred
# For LuxAB
LuxAB_kr1 = 0.28, #inferred
LuxAB_k41 = 0.18, #inferred
LuxAB_k42 = 84.5, #inferred 
LuxAB_k43 = 69.7, #inferred
LuxAB_K_F = 17.1, # Kin. parameter for product inhibition!!!
# For vLuxEC
LuxEC_kr1 = 0.042, # inferred
LuxEC_k61 = 90.92, #inferred
LuxEC_k62 = 95.28, #inferred
LuxEC_k63 = 24.35, #inferred
LuxEC_k64 = 76.51, #inferred
# Degradation constant  (inferred)
Tau = 0.063, # For RCHO and RCOOH
E_LuxAB = 0.32, # uM (estimated from Neils data), divide luxAB V-M by this amount (since that was estimated as V_m*E_LuxAB)
Gam_L  =  0.0063  # /min   0.2
)

ystarts = c(FMNH2=88,RCOOH=115.5,RCHO=115.5,Lux=0)

# Some test values for switch function

time = seq(0,999,1)
PT_time_points = seq(0,1000,1)
PT_heights = c(rep(0,100),rep(1,100),rep(0,801))


# DEFINE MODEL FUNCTION

lux_extended_sub_final = function(time,PT_time_points,PT_heights,params,starts)  {
	
lux3_extended = with(as.list(params),function(t,y,params){
	
vFre = (Vm_Fre*((Conc_F - y["FMNH2"]))*Conc_NADPH) /(Fre_kr1*Fre_k32 + Fre_k32*((Conc_F - y["FMNH2"]))+ Fre_k31*Conc_NADPH + ((Conc_F - y["FMNH2"]))*Conc_NADPH)

vLuxAB = y["Lux"]*((Vm_LuxAB/E_LuxAB)*y["FMNH2"]*Conc_O2*y["RCHO"]) /((LuxAB_kr1*y["RCHO"] + LuxAB_k41*Conc_O2*y["RCHO"] + LuxAB_k43*y["FMNH2"]*Conc_O2 + LuxAB_k42*y["FMNH2"]*y["RCHO"] + y["FMNH2"]*Conc_O2*y["RCHO"])*(LuxAB_K_F + (Conc_F - y["FMNH2"])))
            
   
vLuxEC = (y["Lux"]*Vm_LuxEC)*(y["RCOOH"])*Conc_NADPH*Conc_ATP /(LuxEC_kr1*LuxEC_k62*Conc_NADPH + LuxEC_k61*Conc_ATP*Conc_NADPH + LuxEC_k62*(y["RCOOH"])*Conc_NADPH + LuxEC_k63*(y["RCOOH"])*Conc_ATP*Conc_NADPH + 2*LuxEC_k64*(y["RCOOH"])*Conc_ATP + (y["RCOOH"])*Conc_ATP*Conc_NADPH)

 p0 = Tau*Conc_R
 p1 = 0

# Differential Equations..
dyFMNH2 = vFre - vLuxAB #FMNH2
dyRCOOH =  p1 -1.0*vLuxEC + 1.0*vLuxAB - Tau*y["RCOOH"] # RCOOH   
dyRCHO =  p0 + vLuxEC - vLuxAB - Tau*y["RCHO"]  #RCHO
dyLux =  approx(PT_time_points,PT_heights,t,rule=2)$y  - Gam_L*y["Lux"]   # Generic equation representing lux proteins

 		
list(c(dyFMNH2,dyRCOOH,dyRCHO,dyLux))})
 
# call the ode solver
out = as.data.frame(lsoda(ystarts,time,lux3_extended,params))
 	
# ============================ Output Velocity Calculation ==============

# DJS- this section happens after the ODE solver is called - note the way that it refers to Y, the output of the ODEs.
vLuxAB_out = with(as.list(params),out$Lux*((Vm_LuxAB/E_LuxAB)*out$FMNH2*Conc_O2*out$RCHO) /((LuxAB_kr1*LuxAB_k42*out$RCHO + LuxAB_k41*Conc_O2*out$RCHO + LuxAB_k43*out$FMNH2 *Conc_O2 + LuxAB_k42*out$FMNH2*out$RCHO + out$FMNH2*Conc_O2*out$RCHO)*(LuxAB_K_F + (Conc_F - out$FMNH2 ))))

out$vLuxAB = vLuxAB_out
 	
out
}
 
#======Normalization of data ===========
# All experimental data from one promoter should be normalized by dividing by the maximim light output observed in that promoter



##################
# MCMC Part
##################

MCMC_Prom_Activity = function(T,light,NumIter=10000,frac=0.35,NORMFAC = 10000,Beta=1000,Sigma=0.01,Const_h=1) {

	
# much bigger value of Beta to enforce strong prior on Lux protein abundance	
	
# frac represents the proportion of inferred heights relative to the number of data points


# ===========setting up random stream =========
set.seed(as.numeric(Sys.time()))


# ========Initialize the step positions and heights ===================

# Initializing positions
K = as.integer(frac*length(T))
S_pos = seq(T[1],T[length(T)],length.out=K)
S_ht = rep(0,K)   # Max step heights


# Initializing heights
#height based upon that data,.,,,
for (j in 1:K) {
	S_ht[j] = mean(light)/NORMFAC
} 


# Other Variables used below
h_accepted = 0    
h_tried = 0 
p_accepted = 0     
p_tried = 0


# initialize output matrices.
out_heights=matrix(0,nrow=NumIter,ncol=K)

out_h_acc = rep(0,NumIter)

LikDiff=rep(0,NumIter) # keep these to track relative importance of likelihood and prior
Prop = rep(0,NumIter) # keep these to track relative importance of likelihood and prior

# =========================================================================
# ================== Main RJMCMC loop =====================================
# =========================================================================

# Initial Likelihood calc
 temp = lux_extended_sub_final(T,S_pos,S_ht,luxparams,ystarts)#Call to ODE solver
 SSQ_old = sum((light-temp$vLuxAB)^2)
 Lik = (-0.5/Sigma^2)*SSQ_old #  Calculating summ of sqrs (May be make a function of it)
 
  Lik_h = Lik;
 Lik_p = Lik;


for (r in 1:NumIter) {
	
	if (r %% 1000 == 0) {print(r)}
	        
    
    #================== Standard (M-H) MCMC moves =========================
    # looks like we are doing standard MCMC moves with prob=1 in each
    # iterations..
    
    
    #=========== change a step height =====================================

	K_idx = as.integer(runif(1,max=K))+1 # select the index
	S_ht_new = S_ht[K_idx]*exp(Const_h*(runif(1)-0.5)) # Propose new height
   S_ht_old = S_ht[K_idx]
 	
	

if (K_idx==1) {
	Prop[r] = -S_ht_new*Beta - S_ht[2]/S_ht_new + S_ht_old*Beta + S_ht[2]/S_ht_old	+ log(S_ht_old/S_ht_new)
}
else if (K_idx == K) {
	Prop[r] = -S_ht_new/S_ht[K-1] + S_ht_old/S_ht[K-1]	
}
else {
	Prop[r] = -S_ht_new/S_ht[K_idx-1] - S_ht[K_idx+1]/S_ht_new + S_ht_old/S_ht[K_idx-1] + S_ht[K_idx+1]/S_ht_old	 + log(S_ht_old/S_ht_new)
}

 

   
#if (accept) {# not yet rejected 
    S_ht[K_idx] = S_ht_new # Just for Lik Calculation
    # Calculate (Log) Likelihood
 	temp = lux_extended_sub_final(T,S_pos,S_ht,luxparams,ystarts)#Call to ODE solver
    temp_ssq = sum((light-temp$vLuxAB)^2) # THIS CAN BE IMPROVED BY UPDATING SSQ ONLY FOR THE POINT THAT CHANGES
    Lik = (-0.5/Sigma^2)*temp_ssq #  Calcuating summ of sqrs (May be make a function of it)
    LikDiff[r] = Lik - Lik_h     
    log_Acc_Prob = LikDiff[r] + Prop[r]
	Acc_Prob = exp(max(-30.0,min(0.0,log_Acc_Prob)))
	
    if(runif(1) < Acc_Prob || r==1) {
        Lik_h = Lik
        SSQ_old = temp_ssq
        out_h_acc[r] = 1
    }
    	else {
        S_ht[K_idx] = S_ht_old  
    }
#}
    
    out_heights[r,] = S_ht
   
    
} # end of main loop
list(S_pos,out_heights,out_h_acc,LikDiff,Prop)
} # end of function definition

