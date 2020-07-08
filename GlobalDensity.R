##Global Density
glob <- global_metrics_Majambo_s_MacBook_Pro_Majambo_s_MacBook_Pro
t.glob<-glob[c(6:170)]
t2.glob<-as.data.frame(sapply(t.glob, rescale, to = c(0, 100)))

b <- ggtern(data = glob,
            aes(transi15IF, transi15OF, transi15LR)) + 
  stat_density_tern(geom="polygon",#color='gray',
                    n=500,h=0.1,expand = 1,
                    #bins=50,
                    base='identity',
                    aes(fill   = ..level..),
                    na.rm = TRUE) + 
  #geom_point(color="black",size=1,shape=21) +
  #geom_text(aes(label=id),size=3) + 
  #geom_confidence_tern(contour = TRUE, n = 500, h = NULL, na.rm = FALSE, breaks = c(0.5, 0.9, 0.95),show.legend = NA, 
  #color='black',linetype = "dashed") +
  scale_fill_gradientn(colours = c("blue", "green", "red"),limits=c(0,40),breaks=c(0,10,20,30,40), 
                       labels=c(0,10,20,30,40)) + 
  scale_color_gradientn(colours = c("blue", "green", "red"),limits=c(0,40),breaks=c(0,10,20,30,40), 
                        labels=c(0,10,20,30,40)) + 
  theme_rgbw(base_size = 15) +
  theme(plot.margin=grid::unit(c(0,-5,0,0), "cm")) +
  #theme(legend.justification=c(0,1), legend.position=c(0,1)) +
  #theme_gridsontop() + 
  guides(fill = guide_colorbar(order=1),color="none") + 
  labs(title = "Transient",
       x = "", xarrow = "Indegree", y = "", yarrow = "Outdegree",
       z = "", zarrow = "Local retension", size = 14) +
  theme(plot.title = element_text(hjust = 0.5, vjust= -8, size = 15)) + theme(legend.position = "none")
#labs(title= "Ternary Plot and Filled Contour",fill = "Value, V")


#base + theme(plot.title = element_text(hjust = 0.5, size = 14)) + theme(legend.position = "none")


b1 <- ggtern(data = glob,
             aes(resid15IF, resid15OF, resid15LR)) + 
  stat_density_tern(geom="polygon",#color='gray',
                    n=500,h=0.1,expand = 1,
                    #bins=50,
                    base='identity',
                    aes(fill   = ..level..),
                    na.rm = TRUE) + 
  #geom_point(color="black",size=1,shape=21) +
  #geom_text(aes(label=id),size=3) + 
  #geom_confidence_tern(contour = TRUE, n = 500, h = NULL, na.rm = FALSE, breaks = c(0.5, 0.9, 0.95),show.legend = NA, 
  #color='black',linetype = "dashed") +
  scale_fill_gradientn(colours = c("blue", "green", "red"),limits=c(0,40),breaks=c(0,10,20,30,40), 
                       labels=c(0,10,20,30,40)) + 
  scale_color_gradientn(colours = c("blue", "green", "red"),limits=c(0,40),breaks=c(0,10,20,30,40), 
                        labels=c(0,10,20,30,40)) + 
  theme_rgbw(base_size = 15) +
  theme(plot.margin=grid::unit(c(0,0,0,-4), "cm")) +
  #theme(legend.justification=c(0,1), legend.position=c(0,1)) +
  #theme_gridsontop() + 
  guides(fill = guide_colorbar(order=1),color="none") + 
  labs(title = "Resident",
       x = "", xarrow = "Indegree", y = "", yarrow = "Outdegree",
       z = "", zarrow = "Local retension", size = 14) +
  theme(plot.title = element_text(hjust = 0.5, vjust= -8, size = 15)) + theme(legend.position = "none")
#labs(title= "Ternary Plot and Filled Contour",fill = "Value, V")


#b1 + theme(plot.title = element_text(hjust = 0.5, size = 14)) + theme(legend.position = "none")


b2 <- ggtern(data = glob,
             aes(crypto5IF, crypto5OF, crypto5LR)) + 
  stat_density_tern(geom="polygon",#color='gray',
                    n=500,h=0.1,expand = 1,
                    #bins=50,
                    base='identity',
                    aes(fill   = ..level..),
                    na.rm = TRUE) + 
  #geom_point(color="black",size=1,shape=21) +
  #geom_text(aes(label=id),size=3) + 
  #geom_confidence_tern(contour = TRUE, n = 500, h = NULL, na.rm = FALSE, breaks = c(0.5, 0.9, 0.95),show.legend = NA, 
  #color='black',linetype = "dashed") +
  scale_fill_gradientn(colours = c("blue", "green", "red"),limits=c(0,40),breaks=c(0,10,20,30,40), 
                       labels=c(0,10,20,30,40)) + 
  scale_color_gradientn(colours = c("blue", "green", "red"),limits=c(0,40),breaks=c(0,10,20,30,40), 
                        labels=c(0,10,20,30,40)) + 
  theme_rgbw(base_size = 15) + 
  theme(plot.margin=grid::unit(c(0,-5,0,0), "cm")) +
  #theme(legend.justification=c(0,1), legend.position=c(0,1)) +
  #theme_gridsontop() + 
  guides(fill = guide_colorbar(order=1),color="none") + 
  labs(title = "Cryptobenthic",
       x = "", xarrow = "Indegree", y = "", yarrow = "Outdegree",
       z = "", zarrow = "Local retension", size = 14) + 
  theme(plot.title = element_text(hjust = 0.5, vjust= -8, size = 15)) + theme(legend.position = "none")
#labs(title= "Ternary Plot and Filled Contour",fill = "Value, V")


#base + theme(plot.title = element_text(hjust = 0.5, size = 14)) + theme(legend.position = "none")


b3 <- ggtern(data = glob,
             aes(pare5IF, pare5OF, pare5LR)) + 
  stat_density_tern(geom="polygon",#color='gray',
                    n=500,h=0.1,expand = 1.1,
                    bins=500,
                    base='identity',
                    aes(fill   = ..level..),
                    na.rm = TRUE) + 
  #geom_point(color="black",size=1,shape=21) +
  #geom_text(aes(label=id),size=3) + 
  #geom_confidence_tern(contour = TRUE, n = 500, h = NULL, na.rm = FALSE, breaks = c(0.5, 0.9, 0.95), 
  #color='black',linetype = "dashed") +
  scale_fill_gradientn(colours = c("blue", "green", "red")) + 
  scale_color_gradientn(colours = c("blue", "green", "red")) + 
  theme_rgbw(base_size = 15) + 
  theme(plot.margin=grid::unit(c(0,0,0,-2), "cm")) +
  theme(legend.justification=c(0,1),legend.position = "right") +
  #theme(legend.justification=c(0,1), legend.position=c(0,2)) +
  #theme_gridsontop() + 
  guides(fill = guide_colorbar(order=1),color="none") + 
  labs(title = "Parental",
       x = "", xarrow = "Indegree", y = "", yarrow = "Outdegree",
       z = "", zarrow = "Local retension", size = 14,fill = "Value") +
  theme(plot.title = element_text(hjust = 0.5, vjust= -8, size = 15))  
#labs(title= "Ternary Plot and Filled Contour",fill = "Value, V")

grid.arrange(b, b1, b2, b3, ncol = 2, nrow = 2) 


#dev.copy(pdf,"whatever.pdf")
#dev.off()

#b1 + theme(plot.title = element_text(hjust = 0.5, size = 14)) + theme(legend.position = "none")

