library(stats)
library(fields)
library(akima)
library(viridis)
library(caret)
library(randomForest)
library(ape)
library(gbm)
library(dismo)
library(spdep)
library(MASS)
library(corrplot)
library(geosphere)
library(pdp)

stdize<-function(x,r=1){
  st<-round((x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)),r)
  return(min(x,na.rm=T)+st*(max(x,na.rm=T)-min(x,na.rm=T)))
}


set.seed(04081982)

dir.create('FIGURES')
dir.create('RESULTS')





#####Introduction-cost analysis
costs<-read.csv('cost_table.csv',header=T)
coral.loss.km <- 11700 # sq km lost from 2009-2018
coral.loss.ha <- coral.loss.km * 100

pc.restor.vec <- seq(1,100,1) # percentage restored vector
prop.restor.vec <- pc.restor.vec/100
coral.area.restor <- coral.loss.ha * prop.restor.vec

coral.rest.costpA.lo <- 10000 # USD per hectare
coral.rest.costpA.up <- 5000000 # USD per hectare
med_tot<-median(c(coral.rest.costpA.lo,coral.rest.costpA.up))
costs<-rbind(costs,c('total',0,med_tot,0,coral.rest.costpA.lo,coral.rest.costpA.up))

costs$Median<-as.numeric(costs$Median)
costs$Minimum<-as.numeric(costs$Minimum)
costs$Maximum<-as.numeric(costs$Maximum)


pdf('./FIGURES/fig1.pdf',height=8,width=16)
par(mfrow=c(2,4))
options(scipen=5000)
x<-coral.area.restor/100

for (i in 1:nrow(costs)){
  plot(0,0,type='n',
       ylim=c(min(coral.area.restor)*min(costs$Minimum)/10**9,max(coral.area.restor)*max(costs$Maximum)/10**9),log='y',
       xlim=c(min(x),max(x)),
       ylab='',
       xlab='',las=1,cex.axis=1.2,cex.lab=1.2)
  
  
  axis(side = 3, at = seq(0,max(x),length.out=11), labels = seq(0,100,length.out=11),
       label='aaa') 
  mtext("restored reef (% of reef lost in 2009-2018)", side = 3, line = 2.3,cex=1.2)
  mtext(expression(paste("restored reef (km"^"2", ")")), side = 1, line = 2.3,cex=1.2)
  mtext("restoration costs (US$ billions)", side = 2, line = 3.4,cex=1.2)
  grid(nx = NULL, ny = NULL,
       lty = 2,      # Grid line type
       col = "gray", # Grid line color
       lwd = 0.5)  
  
  text(5000,0.1,costs[i,1],cex=1.5)
  
  coral.restor.cost.lo <- coral.area.restor * costs[i,5]
  coral.restor.cost.up <- coral.area.restor * costs[i,6]
  coral.restor.cost.md <- coral.area.restor * costs[i,3]
  
  
  y_lo<-coral.restor.cost.lo/10**9
  y_up<-coral.restor.cost.up/10**9
  
  
  polygon(c(x,rev(x)),c(y_lo,rev(y_up)),col=viridis(nrow(costs),alpha=0.3)[i],border = NA)
  
  segments(max(x)/10,0.0001,max(x)/10,y_up[which(x==max(x)/10)],lty=2,lwd=0.8)
  segments(max(x)/10,y_up[which(x==max(x)/10)],-2000,y_up[which(x==max(x)/10)],lty=2,lwd=0.8)
  segments(max(x)/10,y_lo[which(x==max(x)/10)],-2000,y_lo[which(x==max(x)/10)],lty=2,lwd=0.8)
  lines(x,coral.restor.cost.md/10**9,col=viridis(nrow(costs))[i],lwd=2)
  
  
}

dev.off()


##############################################
###Analysis 1 - Choice of restored sites

reef<-read.csv('dataset_all_reefs.csv',header=T)
colnames(reef)
data<-reef[,c("restored","lat","lon","remoteness","gravity","cumulative_impacts_mean","cumulative_impacts_trend","baa_n","baa_mean","baa_trend","coral_diversity")]
data<-data[data$gravity>0,]
data<-aggregate(data$restored~.,data=data,FUN='sum')
data<-cbind(data[,ncol(data)],data[,-ncol(data)])


pdf('./FIGURES/ED_fig2.pdf',height=6,width=6)
par(mfrow=c(2,4))
for (i in 4:ncol(data)){
  hist(data[,i],xlab=colnames(data)[i],main='')  
}
dev.off()


M<-cor(data[,4:ncol(data)])**2
diag(M)<-0
pdf('./FIGURES/ED_fig6.pdf',width=4.5,height=4.5)
corrplot(M, method = 'square', order = 'FPC', 
         type = 'lower', diag = FALSE)


dev.off()

####code to identify variables with correlation larger than 0.8
sig_cor<-which(M>0.7,arr.ind=T)
for (i in 1:nrow(sig_cor)){
  print (c(colnames(M)[sig_cor[i,1]],rownames(M)[sig_cor[i,2]]))
}


colnames(data)[1]<-'restored'
data$restored<-1*(data$restored>0)

rest<-which(data$restored==1)
non_rest<-which(data$restored==0)

d <- as.matrix(dist(cbind(data[,2],
                          data[,3])))



var_names<-as.character(colnames(data)[4:ncol(data)])


#############create 1000 sets of points at 100km distance with no autocorrelation
ucoords<-cbind(data$lon,data$lat)
dist<-distm(ucoords, fun = distHaversine)
diag(dist)<-NA
rest_var<-data$restored

tre_d<-150000
val_points<-list()
while (length(val_points)<1000){
    val_p<-1:nrow(dist)
    dist_p<-c()
    while(length(val_p)>0){
      if (length(val_p)==1){
        rand_p<-val_p
      } else {
        rand_p<-sample(val_p,1)}
      dist_p<-c(dist_p,rand_p)
      val_p<-setdiff(val_p,c(which(dist[rand_p,]<tre_d),dist_p))
      #print (length(val_p))
    }
    res<-rest_var[dist_p]
    DT<-data.frame(longitude=ucoords[dist_p,1],
                   latitude=ucoords[dist_p,2])
    
    pts = st_as_sf(DT, coords = c("longitude", "latitude"),
                   crs = 4326, agr = "constant")
    nn = knn2nb(knearneigh(pts,1))
    w = listw2U(nb2listw(nn, style="B", zero.policy = TRUE))
    jct_tot_z<-joincount.multi(as.factor(res), w, zero.policy = TRUE)[4,4]#presences
    if (abs(jct_tot_z)<2){
      val_points[[length(val_points)+1]]<-dist_p
      print (length(val_points))}
    }




loc_frac<-0.1

lrates<-c(0.01,0.001,0.0001)
bag_fractions<-c(0.5,0.7,0.8)
tree_complexities<-(1:5)

w_ratio<-sum(data$restored)/nrow(data)
w_all<-rep(w_ratio,nrow(data))
w_all[which(data$restored==1)]<-1-w_ratio

tot_cal_mod_n<-length(lrates)*length(bag_fractions)*length(tree_complexities)*10
  

brt_cal<-data.frame()
for (rep in 1:10){
  for (l_rate in lrates){
    for (bag_fraction in bag_fractions){
      for (tree_complexity in tree_complexities){
        gbm.fit<-NULL
        while(is.null(gbm.fit)){
          rand_rep<-sample(1:1000,1)
          rest_val<-intersect(rest,val_points[[rand_rep]])
          non_rest_val<-intersect(non_rest,val_points[[rand_rep]])  
          w_ratio<-length(rest_val)/length(val_points[[rand_rep]])
          w_all<-rep(w_ratio,nrow(data))
          w_all[which(data$restored==1)]<-1-w_ratio
          train_ids<-c(sample(rest_val,length(rest_val)*0.8),sample(non_rest_val,length(non_rest_val)*0.8))
          test_ids<-setdiff(val_points[[rand_rep]],train_ids)
          train_data<-data[train_ids,]
          test_data<-data[test_ids,]
          w_sub<-w_all[train_ids]
          gbm.fit <- gbm.step(train_data, gbm.x = 2:ncol(train_data),
                              gbm.y = 1, family="bernoulli", max.trees=10000, tolerance = 0.0001,
                              learning.rate = l_rate, bag.fraction=bag_fraction, 
                              tree.complexity = tree_complexity, silent=T, tolerance.method = "auto",
                              site.weights = w_sub,n.trees=50,step.size=10)
          print ('null')
        }
        pred_data<-predict(gbm.fit,test_data,type='response')
        obs<-test_data$restored
        err<-c()
        for (tre in seq(0.001,0.9,0.001)){
          pred<-1*(pred_data>tre)
          tss_a<-sum((pred>0)*(obs>0))
          tss_b<-sum((pred>0)*(obs==0))
          tss_c<-sum((pred==0)*(obs>0))
          tss_d<-sum((pred==0)*(obs==0))
          #TP / (TP + FN)
          tss_sens<-tss_a/(tss_a+tss_c)
          #TN / (TN + FP)
          tss_spec<-tss_d/(tss_b+tss_d)
          tss<-tss_sens+tss_spec-1
          type_II <- 1-tss_sens #false negative; FN/FN+TP
          type_I<-1-tss_spec #false positive; FP/FP+TN
          err<-rbind(err,as.numeric(c(tre,type_I,type_II,tss)))
        }
        
        
        best_tre<-err[which(err[,4]==max(err[,4])),]
        if (!is.null(dim(best_tre))){
          best_tre<-best_tre[min(which(best_tre[,2]==min(best_tre[,2]))),]
        }      
        brt_cal<-rbind(brt_cal,c(rep,l_rate,bag_fraction,tree_complexity,best_tre,sum(test_data$restored),sum(train_data$restored)))
        if (nrow(brt_cal)==1){colnames(brt_cal)<-c('replicate','l_rate','bag_fraction','tree_complexity',
                                                   'tre','FP','FN','TSS','test_rest','train_rest')}
        print(tail(brt_cal,1))
        print(tot_cal_mod_n-nrow(brt_cal))
      }
    }
  }
}



best_pars<-aggregate(brt_cal$TSS~brt_cal$l_rate+brt_cal$bag_fraction+brt_cal$tree_complexity,FUN='mean')
best_pars<-best_pars[which(best_pars[,4]==max(best_pars[,4]))[1],]

#best_pars<-c(0.001,0.5,3)
l_rate<-best_pars[1]
bag_fraction<-best_pars[2]
tree_complexity<-best_pars[3]



brt_res<-data.frame()
for (rep in 1:1000){
  for (bag in c(0.5)){
    rest_val<-intersect(rest,val_points[[rep]])
    non_rest_val<-intersect(non_rest,val_points[[rep]])  
    w_ratio<-length(rest_val)/length(val_points[[rep]])
    w_all<-rep(w_ratio,nrow(data))
    w_all[which(data$restored==1)]<-1-w_ratio
    gbm.fit<-NULL
    while(is.null(gbm.fit)){
      
      train_ids<-c(sample(rest_val,length(rest_val)*bag),sample(non_rest_val,length(non_rest_val)*bag))
      test_ids<-setdiff(val_points[[rep]],train_ids)
      train_data<-data[train_ids,]
      test_data<-data[test_ids,]
      w_sub<-w_all[train_ids]
      
      
      gbm.fit <- gbm.step(train_data, gbm.x = 2:ncol(train_data),
                          gbm.y = 1, family="bernoulli", max.trees=10000, tolerance = 0.0001,
                          learning.rate = l_rate, bag.fraction=bag_fraction, 
                          tree.complexity = tree_complexity, silent=T, tolerance.method = "auto")
    }
    
    pred_data<-predict(gbm.fit,test_data,type='response')
    obs<-test_data$restored
    err<-c()
    for (tre in seq(0.001,0.9,0.001)){
      pred<-1*(pred_data>tre)
      tss_a<-sum((pred>0)*(obs>0))
      tss_b<-sum((pred>0)*(obs==0))
      tss_c<-sum((pred==0)*(obs>0))
      tss_d<-sum((pred==0)*(obs==0))
      #TP / (TP + FN)
      tss_sens<-tss_a/(tss_a+tss_c)
      #TN / (TN + FP)
      tss_spec<-tss_d/(tss_b+tss_d)
      tss<-tss_sens+tss_spec-1
      type_II <- 1-tss_sens #false negative; FN/FN+TP
      type_I<-1-tss_spec #false positive; FP/FP+TN
      err<-rbind(err,as.numeric(c(tre,type_I,type_II,tss)))
    }
    best_tre<-err[which(err[,4]==max(err[,4])),]
    if (!is.null(dim(best_tre))){
      best_tre<-best_tre[min(which(best_tre[,2]==min(best_tre[,2]))),]
      }      
      brt_res<-rbind(brt_res,c(rep,best_tre,bag))
      if (nrow(brt_res)==1){colnames(brt_res)<-c('replicate','tre','FP','FN','TSS','bag')}
      print(tail(brt_res,1))
    }
}

sink('./RESULTS/cross_validation_restoration_site_choice_summary_05.txt')
print(c('Metric','Mean','SD'))
print(c('FP',mean(brt_res$FP,na.rm=T),sd(brt_res$FP,na.rm=T)))
print(c('FN',mean(brt_res$FN,na.rm=T),sd(brt_res$FN,na.rm=T)))
print(c('TSS',mean(brt_res$TSS,na.rm=T),sd(brt_res$TSS,na.rm=T)))
sink()



write.table(brt_res,'./RESULTS/cross_validation_restoration_site_choice_data.csv',
            row.names=F,col.names=T,sep = ',',quote=F)





###final complete model
train_data<-data[,c(1,4:ncol(data))]
w_ratio<-sum(train_data$restored)/nrow(train_data)
site_weights<-rep(w_ratio,nrow(train_data))
site_weights[which(train_data$restored==1)]<-1-w_ratio



gbm.fit<-NULL
while(is.null(gbm.fit)){
  gbm.fit <- gbm.step(train_data, gbm.x = 2:ncol(train_data),
                      gbm.y = 1, family="bernoulli", max.trees=10000, tolerance = 0.0001,
                      learning.rate = l_rate, bag.fraction=bag_fraction, 
                      tree.complexity = tree_complexity, silent=T, tolerance.method = "auto",
                      site.weights = site_weights)
}


var_names_form<-c("remoteness","gravity","cumulative impacts mean","cumulative impacts trend","severe bleaching alert level events","bleaching alert level mean","bleaching alert level trend","coral diversity")


pred_vals<-predict(gbm.fit,train_data,type='response')
y_min<-min(pred_vals)
y_max<-max(pred_vals)

pdf('./FIGURES/fig2.pdf',width=10,height=5)

par(mfrow=c(2,4))
count<-1
for (var_n in match(gbm.fit$contributions$var,var_names)){

var_name<-var_names[var_n]
pred.grid<-data.frame(seq(min(train_data[,var_n+1]),max(train_data[,var_n+1]),length.out=10000))
colnames(pred.grid)<-var_name

part<-partial(gbm.fit,pred.var=var_name,pred.grid,train=train_data,n.trees=gbm.fit$n.trees,plot=F,type='classification',prob=T)

rel_inf<-round(gbm.fit$contributions$rel.inf[count],1)
count<-count+1
plot(part[,1],part[,2],ylab='probability',
     xlab=paste0(var_names_form[var_n],'\n(',rel_inf,'%)'),
     type='n',
     ylim=c(y_min,y_max),
     las=1,
     cex.axis=1,cex.lab=1
     )
points(train_data[,var_n+1],pred_vals,pch=16,cex=0.5,col="#00007010")
lines(part[,1],part[,2],lwd=1.5,col='darkorange')

rug(train_data[,var_n+1],side=3,col="#00007050",ticksize=-0.1)
}
dev.off()


###2) Determinants of restoration success
#data from literature of coral mortality vs. time
#percentage of yearly mortality for colonies starting from 4-7 cm of size
mort<-c(0,39,39,39,20,20,20,0,0,0,0)
surv<-cumprod((100-mort)/100)*100

t<-seq(0,12*10,12)

surv_data<-data.frame(t=t,surv=surv)
surv_curve <- nls(surv ~ SSasymp(t, yf, y0, log_alpha),data=surv_data)


#simulate restored coral survival associated to post monitoring time (up to 5 years)
rest_surv<-runif(10000,0,100)
t<-data.frame(t=runif(10000,0,12*5)) #


df<-read.csv('dataset_complete.csv',header=T)
df$rest_gen_n[is.na(df$rest_gen_n)]<-1

df<-df[,3:ncol(df)]
df<-df[df$gravity>=0,]

pdf('./FIGURES/ED_fig3.pdf',height=6,width=6)
par(mfrow=c(3,3))
for (i in 10:ncol(df)){
  hist(df[,i],xlab=colnames(df)[i],main='',las=1,cex.lab=1.2,cex.axis=1.2)  
}
dev.off()

pdf('./FIGURES/ED_fig7.pdf',height=5,width=5)
M<-cor(df[,10:ncol(df)])
corrplot(M, method = 'square', order = 'FPC', 
         type = 'lower', diag = FALSE)
dev.off()


t<-data.frame(t=df$Post_monitoring_length)

#derive expected survival from the above relationship
exp_rest_surv<-as.vector(predict(surv_curve,t,data=t))

#measure success as % deviation observed vs expected survival

succ<-100*(exp_rest_surv-df$survival)/exp_rest_surv
succ<-log(1+100-succ)


df$survival<-succ
df<- df[!is.na(df$survival),]
df<- df[!is.na(df$Longitude),]
df<- df[!is.na(df$Latitude),]
df<-df[df$cumulative_impacts_mean>0,]
df$Longitude<-df$Longitude+runif(nrow(df))/100000000


ucoords<-cbind(df$Longitude,df$Latitude)
dist<-distm(ucoords, fun = distHaversine)
diag(dist)<-NA
surv_var<-df$survival


####initial evaluation
res<-surv_var
DT<-data.frame(longitude=ucoords[,1],
               latitude=ucoords[,2])

pts = st_as_sf(DT, coords = c("longitude", "latitude"),
               crs = 4326, agr = "constant")
nn = knn2nb(knearneigh(pts,1))
w = nb2listw(nn, style="W")
moran_p<-moran.test(res, w,alternative = "greater")$p.value


##there are signs of weak correlation; let's take a 50km buffer for safety

tre_d<-1
val_points<-list()
while (length(val_points)<1000){
  val_p<-1:nrow(dist)
  dist_p<-c()
  while(length(val_p)>0){
    if (length(val_p)==1){
      rand_p<-val_p
    } else {
      rand_p<-sample(val_p,1)}
    dist_p<-c(dist_p,rand_p)
    val_p<-setdiff(val_p,c(which(dist[rand_p,]<tre_d),dist_p))
    #print (length(val_p))
  }
  res<-surv_var[dist_p]
  DT<-data.frame(longitude=ucoords[dist_p,1],
                 latitude=ucoords[dist_p,2])
  
  pts = st_as_sf(DT, coords = c("longitude", "latitude"),
                 crs = 4326, agr = "constant")
  nn = knn2nb(knearneigh(pts,1))
  w = nb2listw(nn, style="W")
  moran_p<-moran.test(res, w,alternative = "greater")$p.value
  if (moran_p>0.05){
    val_points[[length(val_points)+1]]<-dist_p
    print (length(val_points))}
}



##first calibration of BRT parameters

var_names<-as.character(colnames(df)[6:ncol(df)])


brt_cal<-c()
for (rep in 1:10){
  for (l_rate in lrates){
    for (bag_fraction in bag_fractions){
      for (tree_complexity in tree_complexities){
        gbm.fit<-NULL
        while (is.null(gbm.fit)){
          rand_rep<-sample(1:1000,1)
          train_ids<-sample(val_points[[rand_rep]],round(length(val_points[[rand_rep]])*0.8))
          test_ids<-setdiff(val_points[[rand_rep]],train_ids)
          train_data<-df[train_ids,c(5:ncol(df))]
          test_data<-df[test_ids,c(5:ncol(df))]
          
          gbm.fit <- gbm.step(train_data, gbm.x = 2:ncol(train_data),
                            gbm.y = 1, family="gaussian", max.trees=10000, tolerance = 0.0001,
                            learning.rate = l_rate, bag.fraction=bag_fraction,
                            tree.complexity = tree_complexity, silent=T, tolerance.method = "auto")
          }
        pred<-predict(gbm.fit,test_data)
        cv_R2<-cor(pred,test_data$survival)**2
        plot(pred,test_data$survival,main=cv_R2)
        brt_cal<-rbind(brt_cal,c(rep,l_rate,bag_fraction,tree_complexity,cv_R2))
        print(tail(brt_cal,1))
          }
        }
       }
    }


best_pars<-aggregate(brt_cal[,5]~brt_cal[,2]+brt_cal[,3]+brt_cal[,4],FUN='mean')
best_pars<-best_pars[which(best_pars[,4]==max(best_pars[,4])),]
l_rate<-best_pars[1]
bag_fraction<-best_pars[2]
tree_complexity<-best_pars[3]


####final cross validation models
gof<-c()
pred_vs_obs_all<-c()
for (rep in 1:1000){
  gbm.fit<-NULL
  while (is.null(gbm.fit)){
    rand_rep<-sample(1:1000,1)
    train_ids<-sample(val_points[[rand_rep]],length(val_points[[rand_rep]])*0.8)
    test_ids<-setdiff(val_points[[rand_rep]],train_ids)
    train_data<-df[train_ids,c(5:ncol(df))]
    test_data<-df[test_ids,c(5:ncol(df))]
    
    gbm.fit <- gbm.step(train_data, gbm.x = 2:ncol(train_data),
                        gbm.y = 1, family="gaussian", max.trees=10000, tolerance = 0.0001,
                        learning.rate = l_rate, bag.fraction=bag_fraction,
                        tree.complexity = tree_complexity, silent=T, tolerance.method = "auto")
help(gbm.step)
  }
  pred<-predict(gbm.fit,test_data)        
  cv_R2<-cor(pred,test_data$survival)**2
  pred_vs_obs_all<-rbind(pred_vs_obs_all,cbind(pred,test_data$survival))      
  gof<-c(gof,cv_R2)
  print(c(rep,tail(gof,1)))
  }


sink('cross_validation_only_restored_summary.txt')
print (c('mean_R2','SD_R2'))
print (c(mean(gof,na.rm=T),sd(gof,na.rm=T)))
sink()


write.table(as.matrix(gof),'cross_validation_only_restored_data.csv',
            row.names=F,col.names=F,sep = ',',quote=F)


pdf('./FIGURES/fig3.pdf',width=4.5,height=4.5)
plot(pred_vs_obs_all[,1],pred_vs_obs_all[,1],
     xlab='predicted',ylab='observed',
     las=1,cex.axis=1.2,cex.lab=1.2,
     pch=16,col='darkblue')
dev.off()


sink('./RESULTS/cross_validation_only_restored_summary.txt')



###check different Alert levels


###
anom<-read.csv('anomalies.csv',header=T)
anom$year_first_anomaly[anom$year_first_anomaly==0]<-2100
anom_res<-anom[anom$restored==1,]
anom_non_res<-anom[anom$restored==0,]

sink('summary_future_dhw20.txt')
mean(anom_res$number_of_anomalies)
sd(anom_res$number_of_anomalies)
mean(anom_res$year_first_anomaly)
sd(anom_res$year_first_anomaly)
sum(anom$number_of_anomalies>0)/nrow(anom)
sum(anom_res$number_of_anomalies>0)/nrow(anom_res)



mean(anom_non_res$number_of_anomalies)
sd(anom_non_res$number_of_anomalies)
mean(anom_non_res$year_first_anomaly)
sd(anom_non_res$year_first_anomaly)
sink()



pdf('time_vs_n_future_bleaching_events.pdf',width=5,height=5)
plot(anom_non_res$number_of_anomalies,anom_non_res$year_first_anomaly,pch=16,
     col='#90909090',cex=0.8,
     xlab=expression("number of years\nwith max(DHW)">=20),
     ylab=expression("first year\nwith max(DHW)">=20),
    cex.lab=1.2,cex.axis=1.2,las=1)


points(anom_res$number_of_anomalies,anom_res$year_first_anomaly,pch=16,
       col='black')

legend(40,2100,c('restored','others'),col=c('black','#90909090'),pch=16,pt.cex=c(1,0.8))
abline(h=mean(anom$year_first_anomaly))
abline(v=mean(anom$number_of_anomalies))
dev.off()


pdf('boxplots.pdf',width=4,height=8)
par(mfrow=c(2,1))
boxplot(log(anom_res$number_of_anomalies),log(anom_non_res$number_of_anomalies),
        names=c('restored','others'),las=1,
        cex.axis=1.2,cex.lab=1.2)

boxplot(anom_res$year_first_anomaly,anom_non_res$year_first_anomaly,
        names=c('restored','others'),las=1,
        cex.axis=1.2,cex.lab=1.2)

dev.off()

mean(anom_res$number_of_anomalies)
mean(anom_non_res$number_of_anomalies)



###per-decade anomalies
library(zoo)
#se<-function(x){return(sd(x)/sqrt(length(x)))}
msd<-function(x,n=10){
  roll_se<-zoo::rollapplyr(x, width = n, FUN = sd, fill = NA,align = 'right')
  return(roll_se)}


y_anom<-read.csv('yearly_anom.csv',header=T)

y<-rollmean(y_anom$anom_n, align='right',10,fill=NA)
y_sd<-msd(y_anom$anom_n)
#
x<-y_anom$year

y<-y[11:length(y)]
y_sd<-y_sd[11:length(y_sd)]
x<-x[11:length(x)]
y_ll<-y-y_sd
y_ul<-y+y_sd

pdf('BE_per_decade.pdf',width=4.5,height=4.5)
plot(x,y,type='n',
     ylim = c(min(y_ll),max(y_ul)),
     xlab='year',
     ylab='mean number of bleaching events per decade',
     las=1,
     cex.axis=1.2,cex.lab=1.2
     
     )

polygon(c(rev(x), x), c(rev(y_ll),(y_ul)),
          col = '#50505050', border = NA)
  
lines(x,y,lwd=1)
dev.off()














plot(y_anom$year,y_anom$anom_n)


