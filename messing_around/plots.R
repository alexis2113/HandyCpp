####
C=B%>% group_by(ST,DATE)%>% summarise(MT=mean(MT)) 
p=ggplot(C, aes(x = DATE, y = , linetype= ST,color=ST, group = ST)) + geom_line(size=1) 
p +scale_linetype_discrete("sole/team",breaks=c("S","T"),labels=c("sole","team")) 
p=last_plot()+ guides(colour=FALSE) 
windowsFonts(TNR = windowsFont(("Times New Roman"))) 
p+theme_few()+theme(text=element_text(family="TNR",size = 14,colour = "black")) 

####
 C=B%>% 
 group_by(model,ST,DATE)%>%  
summarise(Intercept=mean(Intercept)) 
 p=ggplot(C, aes(x = DATE, y = Intercept, linetype= ST,color=ST, group = ST)) + geom_line(size=1)+facet_wrap(~model) p p=p+ylab("Security selection ability")+xlab(NULL) 

p=p+scale_linetype_discrete("sole/team",breaks=c("S","T"),labels=c("sole","team")) 
p=last_plot()+ guides(colour=FALSE) 
 p+theme_few()+theme(text=element_text(family="TNR",size = 12,colour = "black"))
 
 
  ############### 
prow <- plot_grid( p + theme(legend.position="none"),hmp + theme(legend.position="none"), align = 'vh') 
legend_b <- get_legend(p + theme(legend.position="bottom")) 
p2 <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2)) 
p2
