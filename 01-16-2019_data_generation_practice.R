# This is the second practice coding for dissertation simulation study #
# In this coding, only one school level and one student level variable is simualted #
# Note all conditions are estimators are programmed # 
library(BiocParallel)

generate_data <- function(iter, params) {
library(foreign)

n.sim <- 2

datasets = data.frame(replicate = integer(),
                      schid = double(),
                      studentid = double(),
                      y_l = double(),
                      y1_l = double(),
                      y0_l = double(),
                      TE_l = double(),
                      s1_l = double(),
                      s2_l = double(),
                      z2_l = double(),
                      x1_l = double(),
                      v1_l = double(),
                      prob_schl_long = double())
schid = params[iter, "schid"]
studentid = params[iter, "studentid"]
y_l = params[iter, "y_l"]
y1_l = params[iter, "y1_l"]
y0_l = params[iter, "y0_l"]
TE_l = params[iter, "TE_l"]
s1_l = params[iter, "s1_l"]
s2_l = params[iter, "s2_l"]
z2_l = params[iter, "z2_l"]
x1_l = params[iter, "x1_l"]
v1_l = params[iter, "v1_l"]
prob_schl_long = params[iter, "prob_schl_long"]

    for (j in 1:n.sim) {
# Set random seed #
set.seed(1234566 + j)
# Step 1. Generate data #################################################################

# Step 1.1. Generate school level and student level covariates # 
# number of schools = 1, ..., H, number of student per school = 1, ... ,K-h 
# Set the number of schools to be 2000 #
# Generate the number of student per school to be N(200,80), min = 10 and max = 700
# Students per school is based on FL data  
# Generate school-level V1 is continuous N(0, 1)
# Generate student level covariates x1 is continuous N(V1_h, 1)
H=2000
K=as.integer(rnorm(H,200,80))
K[K < 10] = 10
K[K > 700] = 700
v1=rep(NA,H)
v2=rep(NA,H)
x1=NULL
y=NULL

y0=NULL #potential outcome for control units
y1=NULL #potential outcome for treated units
TE=NULL #true treatment effect 
phi0=NULL
phi1=NULL

p2=NULL # school level selection probability
s2=NULL  # s2=1 indicates school selected into sample 
z2=NULL # z2=1 indicates a sampled school was assigned the treatment condition

eta0=eta1=NULL
p1=NULL
s1=NULL

p1_mean=NULL 
s1_mean=NULL
TE_mean=NULL

y1_l=NULL
y0_l=NULL
TE_l=NULL
x1_l=NULL
s1_l=NULL
y_l=NULL
dataset=NULL

  for (h in 1:H){
    # k= # number of student per school, k~N(600, 270), min = 30 and max = 4000
      if (K[h]<10) {
          K[h] = 10
      }
      if (K[h] > 700) {
          K[h] = 700
      }
    v1[h]=rnorm(1,0,1) # v1~N(0,1)
    v2[h]=rbinom(1,1,0.5)
    x1[[h]]=rep(NA,K[h]) # student level variable x1
    x1[[h]]=rnorm(K[h],v1[h],1) # x1~N(v1,1), so v1 is the mean of x1
    # Step 1.2. Simulate potential outcomes for all students in all schools 
    
    y0[[h]]=rep(NA,K[h])
    y1[[h]]=rep(NA,K[h])
    TE[[h]]=rep(NA,K[h])
    
    phi0[h]=pi30+pi31*v1[h] # level 2 intercept as outcome of v1
    phi1[h]=pi40+pi41*v1[h] # level 2 slope as outcome of v1
    
    for (k in 1:K[h]){
      y0[[h]][k]=v1[h]+x1[[h]][k]*(v1[h]) #potential outcome for treatment
      y1[[h]][k]=y0[[h]][k]+phi0[h]+phi1[h]*x1[[h]][k] # potential outcome for control 
      TE[[h]][k]=phi0[h]+phi1[h]*x1[[h]][k] # true individual treatment effect 
    }
    TE_mean[h]=mean(TE[[h]])
    
    
    # Step 1.3.Generate selection probabilities  
    
    # Generate school selection probability and assign treatment groups 
    
    p2[h]=1/(1+exp(-(alpha0+alpha1*v1[h]))) # model for school selection 
    s2[h]=rbinom(1,1,p2[h]) # selection follows a bernoulli distribution 
    if(s2[h]==1) {
      z2[h]=rbinom(1,1,0.5)
    } else {z2[h]=NA}  # assign schools to treatment or control group #
    
    
    # Generate student selection probability 
    # *** This model assumes that every student in the population has their own selection probability, which
    # *** cannot be verified with the data we have 
    
    p1[[h]]=rep(NA,K[h]) # student level probability 
    s1[[h]]=rep(NA,K[h]) # Indicator for student selection
    
    if (s2[h]==1){
      
      eta0[h]=tau00+tau01*v1[h]+v2[h] # level 2 intercept as outcome of v1
      eta1[h]=tau10+tau11*v1[h]+v2[h] # level 2 slope as outcome of v1 
      
      for (k in 1:K[h]){
        p1[[h]][k]=1/(1+exp(-(eta0[h]+eta1[h]*x1[[h]][k]))) # model for student selection 
        s1[[h]][k]=rbinom(1,1,p1[[h]][k]) # student selection follows a bernoulli distribution
      }
      s1_mean[h]=mean(s1[[h]]) # % student selected into sample for school h
    }
    
    # generate outcome variable y for students in the sample schools (s2=1)
    # y=y1 for treatment schools
    # y=y0 for control schools 
    
    y[[h]]=rep(NA,K[h])
    
    if (s2[h]==1 & z2[h]==1) { # if school is in the sample and in the treatment group
      for (k in 1:K[h]){
        if (s1[[h]][k]==1) {
          y[[h]][k]=y1[[h]][k] # and student is sampled, y=y1
        }  else {
          y[[h]][k]= NA # if student is not sampled, y = NA (unobserved)
        }
      }
    }
    else if (s2[h]==1 & z2[h]==0){ # if school is in the sample and in the control group
      for (k in 1:K[h]){
        if (s1[[h]][k]==1) {
          y[[h]][k]=y0[[h]][k]  # and student is sampled, y=y0
        } else {
          y[[h]][k]= NA # if student is not sampled, y = NA (unobserved)
        }
      }
    }
    else {
      for (k in 1:K[h]){
        y[[h]][k]=NA # if school is not sampled, y = NA (unobserved)
      }
    }
  }
  ## mean(s2)  
  ## hist(s1_mean[complete.cases(s1_mean)])
  ## mean(s1_mean[complete.cases(s1_mean)])
  
  # Step 2. Saving all variables into a dataframe ####################################################### 
  
  # Generate school IDs
  schid=rep(seq(1:H),K)
  # Generate student IDs
  studentid=sequence(K)
  
  # flatten all "lists" of variables to make them single variables (this set data to the "long" format) 
  y1_l=unlist(y1)
  y0_l=unlist(y0)
  TE_l=unlist(TE)
  x1_l=unlist(x1)
  s1_l=unlist(s1)
  
  # repeat school-level variables v1,v2,s2,z2 by # students in the school, so that each student has corresponding school level variables  
  v1_l=rep(v1,K)
  s2_l=rep(s2,K)
  z2_l=rep(z2,K)
  
  y_l=unlist(y) # flatten this list and set it to "long" format 
  
  # Combine all variables into one dataset that is ready for estimation 
  # This dataset contains school ID, student ID,
  # outcome for sampled student in sampled schools, 
  # level 1 variables x1 & x2, level 2 variables v1 & v2
  # indicators for school selection and student selection 
  
  schid2=seq(1:H)
  dataset_schl=data.frame(schid2,s2,v1,z2)
  
  # Generate predicte school probability for selection
  fit_schl=glm(s2~v1,data=dataset_schl, family=binomial())
  prob_schl=predict(fit_schl,type="response") # predicted school selection probability 
  prob_schl_long=rep(prob_schl,K) # make it into long format 
  
  # Combine all variables into a single dataset and write to file 
        dataset=data.frame(cbind(replicate=j,schid,studentid,y_l,y1_l,y0_l,TE_l,s1_l,s2_l,z2_l,x1_l,v1_l,prob_schl_long))
        datasets = rbind(datasets, dataset)
        if (j == n.sim) {
            file.name=paste0("dataset",iter,".dta")
            write.dta(datasets,file.name)
        }
    }
return(NULL)
}

# Below are coefficients for the MLM treatment effect # 
# All other coefficients that are in the draft but not listed here are either 0 or 1 #
pi30=1 # unconditional effect of treatment effect (TE) #
pi31=2 # impact of v1 on TE
pi40=2 # impact of x1 on TE
pi41=2 # impact of x1v1 on TE

# Set up school selection parameters. (p is approximately 0.15)
alpha0=-2
alpha1= 0

# Below are coefficiences for the ML logistics regression for student selection 
tau00=-2 # unconditional selection probability 
tau01=0 # impact of v1 on selection
tau10=0 # impact of x1
tau11=0 # impact of x1v1

params <- expand.grid(pi30 = 1,
                      pi31 = c(0.2, 2),
                      pi40 = c(0.2, 2),
                      pi41 = c(0.2, 2),
                      alpha0 = -2,
                      alpha1 = c(0.2, 2),
                      tau00 = -1,
                      tau01 = c(0.2, 2),
                      tau10 = c(0.2, 2),
                      tau11 = c(0.2, 2))

register(MulticoreParam(workers = 12))
#bplapply(1:seq_along(nrow(params)), generate_data, params = params)
generate_data(1, params)
