CairoPNG(file="matirx_scatter.png",width=1200,height=900)
showtext.begin()
ggplot(data=mydata1)+
geom_hline(aes(x=Country,y=Class,yintercept = 1:nrow(mydata1)),size=20,colour="#E4EDF2",alpha=.5)+
geom_vline(aes(x=Country,y=Class,xintercept = 1:nrow(mydata1)),linetype="dashed")+
geom_point(aes(x=Country,y=Class,size=Spend,fill=Spend_fact),shape=21,colour="white")+
scale_fill_manual(values=c("#F9DBD3","#F1B255","#519F46","#41B0C3"))+
scale_size_area(max_size=25)+
scale_x_discrete(position = "top")+
guides(size=FALSE,fill=guide_legend(title="Within category",direction="horizontal"))+
labs(title="How they spend it",subtitle="Househlod spending*,of total,2013 or latest,includes taxes",caption="Source:Eurostat")+
theme_void(base_size=20,base_family="myfont") %+replace%
theme(
      legend.position="top",
      panel.grid.major.x=element_line(linetype="dashed"),      #plot.margin=margin(5,5,5,5,unit="pt"),
      axis.text=element_text(size=15,hjust=0.5),
      plot.title=element_text(size=35,hjust=0,lineheight=1.2),
      plot.caption=element_text(hjust=0,lineheight=1.2)
) 
showtext.end()
dev.off()
