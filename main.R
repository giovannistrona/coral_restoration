library(stats)
library(fields)
library(akima)
library(viridis)
library(caret)
library(randomForest)
library(ape)
library(gbm)
library(dismo)
library(treezy)
library(spdep)
library(MASS)


stdize<-function(x,r=1){
  st<-round((x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)),r)
  return(min(x,na.rm=T)+st*(max(x,na.rm=T)-min(x,na.rm=T)))
}

##############################################
###Analysis 1 - Choice of restored sites

reef<-read.csv('dataset_all_reefs_m.csv',header=T)
colnames(reef)
data<-reef[,-c(11,12)]
rest<-which(data$restored==1)
non_rest<-which(data$restored==0)

d <- as.matrix(dist(cbind(data[,2],
                          data[,3])))



sub_samp_size<-round(length(rest)*0.5)
var_names<-as.character(colnames(data)[4:ncol(data)])

rel_inf<-c()
gof<-c()
p_dep<-c()
sc<-0
while (sc<100){
  rand_sample_rest<-sample(rest,sub_samp_size)
  rand_sample_non_rest<-c()
  for (i in rand_sample_rest){
    rand_sample_non_rest<-c(rand_sample_non_rest,sample(non_rest,1,prob = 1/d[i,non_rest]**2))
  }
  sub_data<-data[c(rand_sample_rest,rand_sample_non_rest),]
  #plot(sub_data[,2],sub_data[,3],pch=19,col=as.factor(sub_data[,1]))
  nn5 = knn2nb(knearneigh(cbind(sub_data[,2]+runif(sub_samp_size*2)/10000,sub_data[,3]),10))
  w = nb2listw(nn5, style="B")
  if (joincount.test(as.factor(sub_data[,1]), w)[[1]]$p.value>0.05){
    brt.fit <- gbm.step(sub_data, gbm.x = 4:ncol(sub_data),
                        gbm.y = 1, family="bernoulli", max.trees=100000, tolerance = 0.0001,
                        learning.rate = 0.001, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
    if (!is.null(brt.fit)){
      summ.fit <- summary(brt.fit)
      rel_inf<-rbind(rel_inf,summ.fit$rel.inf[match(var_names,summ.fit$var)])
      D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
      gof<-c(gof,D2)
      for (v in var_names){
        p_dep<-rbind(p_dep,partial_dependence(brt.fit,var = v))
      }
      sc<-sc+1
      print(c(nrow(rel_inf),D2))
    }
  } else {
    print ('spatial autocorrelation!')}
}




colnames(rel_inf)<-var_names

write.table(rel_inf,'relative_influence_rest_loc.csv',
            row.names=F,col.names=T,sep = ',',quote=F)

write.table(p_dep,'partial_dep_rest_loc.csv',
            row.names=F,col.names=T,sep = ',',quote=F)

write.table(as.matrix(gof),'gof_rest_loc.csv',
            row.names=F,col.names=F,sep = ',',quote=F)



gof<-read.csv('gof_rest_loc.csv',header=F)
head(gof)
mean(gof[,1])
sd(gof[,1])

rel_inf<-read.csv('relative_influence_rest_loc.csv',header=F)
head(rel_inf)
dim(rel_inf)



pdf('brt_rest_sites_var_imp.pdf')
# variable relative importance
rel_inf<-rel_inf[,order(colMeans(rel_inf))]
boxplot(rel_inf,horizontal=T,las=1,outline=F,cex.axis=0.5)

#weigh by D2
rel_inf_w<-t(t(rel_inf)*(gof/100))
rel_inf_w<-rel_inf_w[,order(colMeans(rel_inf_w))]
boxplot(rel_inf_w,horizontal=T,las=1,outline=F,cex.axis=0.5)
dev.off()


top_vars<-rev(colnames(rel_inf_w))
pdf("partial_dependence_rest_loc.pdf",width=15,height=8)
par(mfrow=c(2,4))
for (var in top_vars){
  p_dep_v<-p_dep[p_dep$variable==var,]
  p_dep_v$value<-stdize(p_dep_v$value,2)
  xy<-aggregate(p_dep_v$fitted_function~p_dep_v$value,
                FUN='mean',na.rm=TRUE, na.action=NULL)
  x<-xy[,1]
  y<-xy[,2]
  xy_sd<-aggregate(p_dep_v$fitted_function~p_dep_v$value,FUN='sd',na.rm=T)[,2]
  xy_sd[which(is.na(xy_sd))]<-0
  xy_n<-aggregate(p_dep_v$fitted_function~p_dep_v$value,FUN='length')[,2]
  ci<-qnorm(0.99)*xy_sd/sqrt(xy_n)
  y_min<-min(y-ci)
  y_max<-max(y+ci)
  x_min<-min(x)
  x_max<-max(x)
  plot(x,y,type='n',las=1,xlab=var,ylab='restoration site selection p',
       xlim=c(x_min,x_max),ylim=c(y_min,y_max))
  polygon(c(rev(x), x), c(rev(y-ci),(y+ci)),
          col = plasma(1,alpha=0.5)[1], border = NA)
  lines(x,y,lwd=1.5,col=plasma(1)[1])
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


df<-read.csv('dataset_only_restored_m.csv',header=T)
df<- df[!is.na(df$survival),]

t<-data.frame(t=df$Post_monitoring_length)

#derive expected survival from the above relationship
exp_rest_surv<-as.vector(predict(surv_curve,t,data=t))

#measure success as % deviation observed vs expected survival

succ<-100*(exp_rest_surv-df$survival)/exp_rest_surv
#succ[succ<0]<-0
succ<-log(1+100-succ)
hist(succ)
plot(succ,df$survival)


df$survival<-succ
df<- df[!is.na(df$survival),]
df<- df[!is.na(df$Longitude),]
df<- df[!is.na(df$Latitude),]
df<-df[df$cumulative_impacts_mean>0,]
df$Longitude<-df$Longitude+runif(nrow(df))/100000000


d <- as.matrix(dist(cbind(df$Longitude,df$Latitude)))
d.inv <- 1/d
diag(d.inv) <- 0
Moran.I(df$survival, d.inv)$p.value


####BRT to look for potential model predicting success of restoration based on multiple variables
#colnames(df)
#rf$Transplantation<-as.factor(rf$Transplantation)
#rf$Coral_gardening<-as.factor(rf$Coral_gardening)
# rf$Substrate_enhancement<-as.factor(rf$Substrate_enhancement)


sub_samp_size<-round(nrow(df)*0.85)
var_names<-as.character(colnames(df)[6:ncol(df)])

rel_inf<-c()
p_dep<-c()
gof<-c()
sc<-0
while (sc<100){
  sub_data<-df[sample(1:nrow(df),sub_samp_size),]
  d <- as.matrix(dist(cbind(sub_data[,2],
                            sub_data[,3])))
  d.inv <- 1/d
  diag(d.inv) <- 0
  if (Moran.I(sub_data[,1], d.inv)$p.value>0.05){
    brt.fit <- gbm.step(sub_data, gbm.x = 6:ncol(sub_data),
                      gbm.y = 5, family="gaussian", max.trees=100000, tolerance = 0.0001,
                      learning.rate = 0.001, bag.fraction=0.75, tree.complexity = 2, silent=T, tolerance.method = "auto")
    if (!is.null(brt.fit)){
      summ.fit <- summary(brt.fit)
      rel_inf<-rbind(rel_inf,summ.fit$rel.inf[match(var_names,summ.fit$var)])
      D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
      gof<-c(gof,D2)
      for (v in var_names){
        p_dep<-rbind(p_dep,partial_dependence(brt.fit,var = v))
      }
      sc<-sc+1
      print(c(nrow(rel_inf),D2))
    }
    else {print ('spatial autocor')}
  }
}



colnames(rel_inf)<-var_names
write.table(rel_inf,'relative_influence_success.csv',
            row.names=F,col.names=T,sep = ',',quote=F)

write.table(p_dep,'partial_dep_success.csv',
            row.names=F,col.names=T,sep = ',',quote=F)

write.table(as.matrix(gof),'gof_success.csv',
            row.names=F,col.names=F,sep = ',',quote=F)


gof<-read.csv('gof_success.csv',header=F)
summary(gof[,1])

170/299

pdf('brt_succ_var_imp.pdf')
# variable relative importance
rel_inf<-rel_inf[,order(colMeans(rel_inf))]
boxplot(rel_inf,horizontal=T,las=1,outline=F,cex.axis=0.5)
#weigh by D2

rel_inf_w<-t(t(rel_inf)*(gof/100))
rel_inf_w<-rel_inf_w[,order(colMeans(rel_inf_w))]
boxplot(rel_inf_w,horizontal=T,las=1,outline=F,cex.axis=0.5)
dev.off()



top_vars<-rev(colnames(rel_inf))[1:9]
pdf("partial_dependence_success.pdf",width=10,height=10)
par(mfrow=c(3,3))
for (var in top_vars){
  p_dep_v<-p_dep[p_dep$variable==var,]
  p_dep_v$value<-stdize(p_dep_v$value,2)
  xy<-aggregate(p_dep_v$fitted_function~p_dep_v$value,
                FUN='mean',na.rm=TRUE, na.action=NULL)
  x<-xy[,1]
  y<-xy[,2]
  xy_sd<-aggregate(p_dep_v$fitted_function~p_dep_v$value,FUN='sd',na.rm=TRUE, na.action=NULL)[,2]
  xy_sd[which(is.na(xy_sd))]<-0
  xy_n<-aggregate(p_dep_v$fitted_function~p_dep_v$value,FUN='length')[,2]
  ci<-qnorm(0.99)*xy_sd/sqrt(xy_n)
  y_min<-min(y-ci)
  y_max<-max(y+ci)
  x_min<-min(x)
  x_max<-max(x)
  plot(x,y,type='n',las=1,xlab=var,ylab='restoration success',
       xlim=c(x_min,x_max),ylim=c(y_min,y_max))
  polygon(c(rev(x), x), c(rev(y-ci),(y+ci)),
        col = plasma(1,alpha=0.5)[1], border = NA)
lines(x,y,lwd=1.5,col=plasma(1)[1])
}
dev.off()



###check different Alert levels
a<-read.csv("BAA.csv",header=T)

a$Year<-a$Year-a$Year%%5
mean_baa<-aggregate(1*(a$baa3_n_post>0)~a$Year,FUN='sum')
n_rest<-aggregate(a$baa3_n_post~a$Year,FUN='length')

cairo_pdf('baa_trend.pdf')
plot(mean_baa[,1],100*mean_baa[,2]/n_rest[,2],type='l',
     xlab='year',ylab='% restored sites\nexposed to bleaching alert level \u2265 1',
     las=1,cex.axis=1.5,cex.lab=1.5)
points(mean_baa[,1],100*mean_baa[,2]/n_rest[,2],pch=16,cex=4,col='lightgrey')
text(mean_baa[,1],100*mean_baa[,2]/n_rest[,2],n_rest[,2])
dev.off()






mean_baa<-aggregate(1*(a$baa3_n_post>0)~round(a$Year/5),FUN='sum')
n_rest<-aggregate(a$baa3_n_post~round(a$Year/5),FUN='length')

cairo_pdf('baa_trend.pdf')
plot(mean_baa[,1]*5,100*mean_baa[,2]/n_rest[,2],type='l',
     xlab='year',ylab='% restored sites\nexposed to bleaching alert level \u2265 1',
     las=1,cex.axis=1.5,cex.lab=1.5)
points(mean_baa[,1]*5,100*mean_baa[,2]/n_rest[,2],pch=16,cex=4,col='lightgrey')
text(mean_baa[,1]*5,100*mean_baa[,2]/n_rest[,2],n_rest[,2])
dev.off()



sum(a$baa3_n_post>0)/nrow(a)
summary(a$baa3_n_post[which(a$baa3_n_post>0)])
sd(a$baa3_n_post[which(a$baa3_n_post>0)])
summary(a$baa3_n_post)
sd(a$baa3_n_post)


sum(a$baa4_n_post>0)/nrow(a)
summary(a$baa4_n_post[which(a$baa4_n_post>0)])
sd(a$baa4_n_post[which(a$baa4_n_post>0)])
summary(a$baa4_n_post)
sd(a$baa4_n_post)

