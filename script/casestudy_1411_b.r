 

data = annotation%>%arrange(start1)%>%filter(start2 > 36500000)
ggplot() + 
 
  geom_rect(data=data[-which(duplicated(data$start1)),], mapping=aes( 
    xmin = start1,
    xmax = end1,
    ymin = 0.87,
    ymax =  1.08,
    color = region, fill=region) , alpha=0.5)  +
  geom_rect(data=data[-which(duplicated(data$start1)),],aes(xmin = start2,
                             xmax = end2,
                             ymin = 0.87,
                             ymax =  1.08),
               alpha = 0.2, 
               color = "green",
            fill= "green")+
  geom_segment(data=data,aes( x = start1,
                    xend = (start1+start2)/2,
                    y = 0.87 , 
                    yend = 0.835 )
  )+
  geom_segment(data=data,aes( x = (start1+start2)/2 ,
                    xend = start2, 
                    y = 0.835 , 
                    yend = 0.87 )
  )+
  geom_segment(data = tf_candidates%>%
                 filter() ,
               aes(x = X2,
                   xend = X3,
                   y = 1.1-X7/30,
                   yend = 1.1-X7/30),# bars
               size = 2
  )+
  geom_segment( aes(x = start,
                    xend = end, 
                    y = (nn-strand/100),
                    yend =(nn-strand/100)
  ) ,
  arrow = arrow(length = unit(0.5, "cm")),
  data = gene_plot # arrow gene
  ) +
  geom_text(aes(x = 35985000, 
                y = 1.1-X7/30,
                label = X4 ),
            size = 5,check_overlap = TRUE,
            data =tf_candidates%>%filter() # names of the TF
  ) +
  geom_text(data = gene_plot,
            aes(x = (x_label+30000),
                y = 0.885, 
                label = gene_name,
                vjust=-1),#name gene
            size = 5
  ) +
  scale_color_manual("region",
                     values = c("orange",
                                "red",
                                "yellow",
                                "purple") 
  )+
  scale_fill_manual("region",
                     values = c("orange",
                                "red",
                                "yellow",
                                "purple") 
  )
