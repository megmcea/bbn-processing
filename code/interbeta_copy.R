
#This is a copy of code made by Mike Colvin, USGS
#last updated 8/18/2026

# I have adapted it to work better as a functionalized code that can be applied to any of the sheets I give it.

#load libraries
library(tidyverse)

#read in all expert responses to the bbn structure:
response_files<-list.files(path="data/",pattern="*.csv",full.names=T)
responses_list<-lapply(response_files,read.csv)
names(responses_list)<-str_extract(response_files, "(?<=_)[^_]+(?=\\.csv$)")



#recode map for nodes to make them more manageable:
recode_map_nodes<-c(
  "Combined Costs and Benefits (Favorable / Unfavorable) "="costs_benefits_node",
  "Deployable in a Reasonable Timeframe "="reasonable_deployment",
  "Will the Control Tool Be Successful? "="successful_tool",
  "Social Acceptance "="social_acceptance"
  
)

#process the responses to prepare for generating beta distributions for the interbeta nodes:
response<-responses_list[[1]] #for testing

process_response_bestworst<-function(response){
  
  #grab out just the best and worst panels (these are for the beta distributions)
  response %>% filter(panel=="Best/Worst Case Panels")->best_worst
  
  #split response columns 
  best_worst %>% 
    separate_wider_delim(cols = item, delim  = "-", names=c("node","best_or_worst","type"))->best_worst
  
  #recode for easier working with them
  best_worst<-best_worst %>% 
    mutate(node= recode(node, !!! recode_map_nodes ))
  
  return(best_worst)
  
}

processed_bestworst_responses<-map(responses_list, ~process_response_bestworst(.x))

#Ok now have a list of the best and worst responses processed. Need to apply the function to get beta params

#Bind into single df for later processing:
bound_data_bestworst<-bind_rows(processed_bestworst_responses,.id="respondent")
#==================================================================
# Helper functions -- left these untouched from Mike's code 
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


#now from each expert elicitation, for each node, I need to extract these values:
# 3 states (high, moderate, low) for child node
#these have to come from the expert elicitation forms
best_case<-c(0.95,0.04,0.01)
worst_case<- c(0.05,0.15,0.8)

#making this a functionalized version of mike's code that takes as its input each of the questions posed by each of the experts:
build_cpts<-function(x){ #takes as its argument the node of interest (i.e. social_acceptance)
  
  #uses the data stored in bound_data_bestworst df: 
  bound_data_bestworst %>% 
    filter(node==x) ->node_data  #filter to node of interest
  
  #now have to pull out best and worst values for each expert:
  for (i in unique(node_data$respondent)){
    respondent_dat<-node_data[respondent==i,]
    
  }
  
####  ----- what follows is mike's code, unadulterated: 
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


