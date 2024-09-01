setwd("D:\\Liu_analysis")
REE5_all = read.csv("E:/allmetal_type_clin.csv",row.names = 1)

###BKMR model
REE5_all <- subset(REE5_all, REE5_all$pts == "LC")
covar <- data.matrix(REE5_all[, c("age")])
expos <- data.matrix(REE5_all[, c("logLaNM", "logCeNM", "logPrNM","logNdNM", "logYNM")])
Y <- as.numeric(REE5_all$snv_num+REE5_all$indel_num)

###整体暴露
set.seed(906)
fitpr <- kmbayes (y = Y, Z = expos, X = covar, iter = 1000,
                  verbose = FALSE, varsel = TRUE, 
                  est.h = TRUE,family = "gaussian",
                  control.params = list(r.jump2 = 0.5))

fitpr
TracePlot(fit = fitpr, par = "beta")
TracePlot(fit = fitpr, par = "sigsq.eps")
TracePlot(fit = fitpr, par = "r", comp = 7)

ExtractPIPs(fitpr)


risks.overall <- OverallRiskSummaries(fit = fitpr, qs = seq(0.25, 0.75, by = 0.05), 
                                      q.fixed = 0.5)
risks.overall

p = ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, 
                              ymax = est + 1.96*sd))
r <- p+geom_hline(yintercept = 0, lty = 3, col = "brown") +
  geom_pointrange()+theme_bw()+theme(panel.grid = element_blank())+ylab("Overall effect")+theme(
    axis.text.x = element_text(face="bold",color="black",size = 15),
    axis.text.y = element_text(face="bold",color="black",size = 15),
    axis.title.y = element_text(face="bold",size =18 ),
    axis.title.x = element_text(face="bold",size = 18))+
  scale_size_manual(values = c(100))
r
ggsave(paste0("REE/pictures/BKMR_inframedel_1.pdf"),plot = r, width = 8, height = 5)


##单个因素影响
risks.singvar <- SingVarRiskSummaries(
  fit = fitpr, qs.diff = c(0.25, 0.75), 
  q.fixed = c(0.25, 0.50, 0.75))

q <- ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd, 
                               col = q.fixed, shape = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + coord_flip()+theme_bw()+
  theme(panel.grid = element_blank())+ylab("Single-exposure effect")+theme(
    axis.text.x = element_text(face="bold", color="black",size = 15),
    axis.text.y = element_text(face="bold",color="black", size = 15),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face="bold",size = 18))+
  geom_hline(yintercept = 0, lty = 2, col = "black")+theme(legend.position = c(0.95,0.9),
                                                           legend.justification =c(1,1))+
  scale_color_manual(values = c("#BB191A","#FE8002","#1764AA"))+
  scale_size_manual(values = c(5,10,15))+
  scale_shape_manual(values = c(17, 19, 15))

q

ggsave(paste0("REE/pictures/BKMR_F_indel_2.pdf"),plot = q, width = 8, height = 5)


######

pred.resp.univar <- PredictorResponseUnivar(fit = fitpr)

r <- ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(exposure)")+xlab("exposure")+theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values = c("darkred","grey","grey"))+
  theme(axis.text.x = element_text(face="bold", color="black", size=15))+
  theme(axis.text.y = element_text(face="bold", color="black", size=15))+
  theme(axis.title.x = element_text(face="bold", color="black", size=18))+
  theme(axis.title.y = element_text(face="bold", color="black", size=18))
r
ggsave(paste0("REE/pictures/BKMR_F_indel_3.pdf"),plot = r, width = 8, height = 5)


