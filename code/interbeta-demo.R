
#==================================================================
# Helper functions 
#==================================================================

# Function to calculate Kullback–Leibler divergence to use as an
# objective function to minimize in optimization
kld<-function(P=NULL, Q=NULL) {
    return(sum(P*log(P/Q)))
}

# Calculate multinomial probabilities from beta function
# given an alpha and beta
kld_beta<-function(p=NULL, 
    elicited_values=NULL) {
    alpha<-p[1]
    beta<-p[2]
    # Equal intervals for states
    pr<-seq(0,1,length.out=length(elicited_values)+1)
    # Calculate Kullback–Leibler divergence
    D<-kld(P=elicited_values, Q=diff(pbeta(pr, alpha, beta)))
    return(D)
}



#==================================================================
# Fit beta distributions to best and worst case scenarios 
#==================================================================

# 3 states (high, moderate, low) for child node
best_case<-c(0.95,0.04,0.01)
worst_case<- c(0.05,0.15,0.8)

grd<-expand.grid(alpha=seq(0,15, by=0.1), beta=seq(0,15,by=0.1))
# Grid search to find some reasonable starts to fit the beta 
grd$best<-sapply(seq.int(nrow(grd)), function(x,case_dat) {
    kld_beta(p=c(grd$alpha[x],grd$beta[x]),
        elicited_values=case_dat)
},case_dat=best_case)
# Estimate beta distribution for best case scenario
est_best<- optim(par= grd[which.min(grd$best),c(1:2)], 
    fn=kld_beta,
    elicited_values=best_case,
    method = "BFGS")
# Grid search to find some reasonable starts to fit the beta 
grd$worst<-sapply(seq.int(nrow(grd)), function(x,case_dat) {
    kld_beta(p=c(grd$alpha[x],grd$beta[x]),
        elicited_values=case_dat)
},case_dat=worst_case)
# Estimate beta distribution for worst case scenario
est_worst<- optim(par=grd[which.min(grd$worst),c(1:2)], 
    fn=kld_beta,
    elicited_values=worst_case,
    method = "BFGS")


#==================================================================
# Build cpt table to interpolate probabilities to from
# best and worst case scenarios
#==================================================================

# States for parents node in order
parent_states<-list(
    p1=c("high","moderate","low"),
    p2=c("high","moderate","low"),
    p3=c("high","low"),
    p4=c("high","moderate","low","very low"))
# States for child node in order
child_states<-c("high","moderate","low")


# Weights for 4 parent nodes
weights<-c(p1=0.3,p2=0.2,p3=0.4,p4=0.1)
# Make cpt--state order is important, expand.grid respects order 
# when applying factor levels
cpt<-do.call(expand.grid,parent_states)
# Create a table of row scores for the cpt
row_scores<- cpt
for(i in 1:ncol(row_scores)) {
    # Convert factor level to ordinal scale
    row_scores[,i]<-as.numeric(row_scores[,i])
    y<-row_scores[,i]
    # Normalize y such that the desired state has the 
    # highest score and the least desired state is 0. 
    y<-(max(y)-y)/(max(y)-min(y)) 
    # Weight parent specific row score
    row_scores[,i]<-weights[i]*y
}
# Calculate overall row scores
row_scores$score<-apply(row_scores,1,mean)
# Normalize score for interpolation
row_scores$interp_score<-(row_scores$score-min(row_scores$score))/(max(row_scores$score)-min(row_scores$score))
# Interpolate alpha
row_scores$alpha<-row_scores$interp_score*est_best$par[1]+(1-row_scores$interp_score)*est_worst$par[1]
# Interpolate beta
row_scores$beta<-row_scores$interp_score*est_best$par[2]+(1-row_scores$interp_score)*est_worst$par[2]



#==================================================================
# Calculate the probability for the child node outcomes 
#==================================================================

pr_outcomes<-as.data.frame(matrix(0, nrow=nrow(cpt), ncol=length(child_states)))
names(pr_outcomes)<-child_states

for(i in 1:nrow(row_scores)) {
    # Equal intervals for each state
    pp<-seq(0,1,length.out=length(child_states)+1)
    # Probability for each state
    y<-as.data.frame(matrix(diff(pbeta(pp, row_scores$alpha[i], 
        row_scores$beta[i])),nrow=1))
    pr_outcomes[i,]<-y
}

cpt<-cbind(cpt, pr_outcomes) 
# Parameterized cpt
cpt


