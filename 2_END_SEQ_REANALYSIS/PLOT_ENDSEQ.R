library(rtracklayer)
library(ggplot2)
library(broom)
bg_for <- as.data.frame(import.bw('END_WTC_STAT_for.bigwig'))[,c(2,6)]
bg_rev <- as.data.frame(import.bw('END_WTC_STAT_rev.bigwig'))[,c(2,6)]

colnames(bg_for) <- c("start", "value")
colnames(bg_rev) <- c("start", "value")

oriC <- 3926000
interval <- 1000
adjustment <- oriC - median(bg_for$start)
if (adjustment < 0) {
  bg_for$adjusted_start <- ifelse(bg_for$start > (max(bg_for$start) + adjustment), 
                                  bg_for$start - (max(bg_for$start) + interval), bg_for$start)} else {
                                    bg_for$adjusted_start <- ifelse(bg_for$start < (min(bg_for$start) + adjustment), 
                                                                    bg_for$start + (max(bg_for$start) + interval), bg_for$start)}
if (adjustment < 0) {
  bg_rev$adjusted_start <- ifelse(bg_rev$start > (max(bg_rev$start) + adjustment), 
                                  bg_rev$start - (max(bg_rev$start) + interval), bg_rev$start)} else {
                                    bg_rev$adjusted_start <- ifelse(bg_rev$start < (min(bg_rev$start) + adjustment), 
                                                                    bg_rev$start + (max(bg_rev$start) + interval), bg_rev$start)}

fit_for <- loess(value ~ adjusted_start, data = bg_for, span = 0.1)
fit_rev <- loess(value ~ adjusted_start, data = bg_rev, span = 0.1)
pred_for <- augment(fit_for, newdata = data.frame(adjusted_start = bg_for$adjusted_start))
pred_rev <- augment(fit_rev, newdata = data.frame(adjusted_start = bg_for$adjusted_start))
colnames(pred_for) <- c('start','value')
colnames(pred_rev) <- c('start','value')
rm(bg_for,bg_rev,fit_for,fit_rev,adjustment,interval)

p <- ggplot()+
  geom_line(data=pred_for,aes(x=start, y=value),color="#9eb157",linewidth=0.5)+
  geom_ribbon(data = pred_for, aes(x = start, ymin = 0, ymax = value), fill = "#9eb157", alpha = 0.3) +
  geom_line(data=pred_rev,aes(x=start, y=value),color="#bda7e3",linewidth=0.5)+
  geom_ribbon(data = pred_rev, aes(x = start, ymin = 0, ymax = value), fill = "#bda7e3", alpha = 0.3) +
  geom_vline(xintercept=oriC, linewidth=1,linetype="dashed",color='#414042')+
  scale_y_continuous(limits=c(0,510),
                     expand = expansion(mult = c(0,0)))+
  scale_x_continuous(breaks=c(oriC-2e6,oriC-1e6,oriC,oriC+1e6,oriC+2e6),
                     labels=c('-2MB','-1MB','oriC','+1MB','+2MB'),
                     expand = expansion(mult = c(0,0)))+
  theme_classic()+
  theme(legend.position = 'none', 
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_line(color = '#414042'),
        axis.line = element_line(color = '#414042'),
        axis.text = element_text(color = '#58595b'),
        strip.text = element_blank())
print(p)
pdf(paste0(wd,"stat_ENDSEQ.pdf"),width=2.1,height=1.5,colormodel='rgb',paper = 'A4')
print(p)
dev.off()





bg_for <- as.data.frame(import.bw('END_WTC_STAT_for.bigwig'))[,c(2,6)]
bg_rev <- as.data.frame(import.bw('END_WTC_STAT_rev.bigwig'))[,c(2,6)]

colnames(bg_for) <- c("start", "value")
colnames(bg_rev) <- c("start", "value")

oriC <- 3926000
interval <- 1000
adjustment <- oriC - median(bg_for$start)
if (adjustment < 0) {
  bg_for$adjusted_start <- ifelse(bg_for$start > (max(bg_for$start) + adjustment), 
                                  bg_for$start - (max(bg_for$start) + interval), bg_for$start)} else {
                                    bg_for$adjusted_start <- ifelse(bg_for$start < (min(bg_for$start) + adjustment), 
                                                                    bg_for$start + (max(bg_for$start) + interval), bg_for$start)}
if (adjustment < 0) {
  bg_rev$adjusted_start <- ifelse(bg_rev$start > (max(bg_rev$start) + adjustment), 
                                  bg_rev$start - (max(bg_rev$start) + interval), bg_rev$start)} else {
                                    bg_rev$adjusted_start <- ifelse(bg_rev$start < (min(bg_rev$start) + adjustment), 
                                                                    bg_rev$start + (max(bg_rev$start) + interval), bg_rev$start)}
lf <- sum(subset(bg_for, adjusted_start < oriC)$value)
rf <- sum(subset(bg_for, adjusted_start > oriC)$value)
lr <- sum(subset(bg_rev, adjusted_start < oriC)$value)
rr <- sum(subset(bg_rev, adjusted_start > oriC)$value)
all <- lf + rf + lr + rr
lf_frac <- lf/all
rf_frac <- rf/all
lr_frac <- lr/all
rr_frac <- rr/all
