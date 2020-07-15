#Majambo
#Plot Global scatter on triangle

library(ggplot2)
library(ggtern)
library(tidyverse)
library(patchwork)
library(scales)

#________________________________________________________________________
#Read Data
#________________________________________________________________________

glob <- read.csv('_data/nanscaledTriangle.csv')
t.glob<-glob[c(6:170)]
t2.glob<-as.data.frame(sapply(t.glob, rescale, to = c(0, 100)))

b <- ggtern(data = t2.glob,
            aes(transi15in, transi15out, transi15LR)) + 
  geom_point(color="gray10",size=2,shape=21) +
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


b1 <- ggtern(data = t2.glob,
             aes(resid15in, resid15out, resid15LR)) + 
  geom_point(color="gray10",size=2,shape=21) +
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


b2 <- ggtern(data = t2.glob,
             aes(crypto5in, crypto5out, crypto5LR)) + 
  geom_point(color="gray10",size=2,shape=21) +
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
  geom_point(color="gray10",size=2,shape=21) +
  theme_rgbw(base_size = 15) + 
  theme(plot.margin=grid::unit(c(0,0,0,-4), "cm")) +
  #theme(legend.justification=c(0,1), legend.position=c(0,1)) +
  #theme_gridsontop() + 
  #geom_confidence_tern(contour = TRUE, h = NULL, na.rm = FALSE, breaks = c(0.5, 0.9, 0.95), 
  #                    color='orange',linetype = "dashed") +
  guides(fill = guide_colorbar(order=1),color="none") + 
  labs(title = "Parental",
       x = "", xarrow = "Indegree", y = "", yarrow = "Outdegree",
       z = "", zarrow = "Local retension", size = 14) +
  theme(plot.title = element_text(hjust = 0.5, vjust= -8, size = 15)) + theme(legend.position = "none")
#labs(title= "Ternary Plot and Filled Contour",fill = "Value, V")

grid.arrange(b, b1, b2, b3, ncol = 2, nrow = 2) 


#________________________________________________________________________
#Relative Proportions
#________________________________________________________________________

glob <- read.csv('_data/relative_proportions.csv')

b <- ggtern(data = glob,
            aes(transi15in, transi15out, transi15LR)) + 
  geom_point(color="gray10",size=1,shape=21) +
  theme_rgbw(base_size = 15) +
  theme(plot.margin=grid::unit(c(0,-5,0,0), "cm")) +
  geom_confidence_tern(contour = TRUE,h = NULL, na.rm = FALSE, breaks = c(0.2, 0.5, 0.75),show.legend = NA, 
                       color='orange',size=0.5) +
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
             aes(resid15in, resid15out, resid15LR)) + 
  geom_point(color="gray10",size=1,shape=21) +
  theme_rgbw(base_size = 15) +
  theme(plot.margin=grid::unit(c(0,0,0,-4), "cm")) +
  geom_confidence_tern(contour = TRUE,h = NULL, na.rm = FALSE, breaks = c(0.2, 0.5, 0.75),show.legend = NA, 
                       color='orange',size=0.5) +
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
             aes(crypto5in, crypto5out, crypto5LR)) + 
  geom_point(color="gray10",size=1,shape=21) +
  theme_rgbw(base_size = 15) + 
  theme(plot.margin=grid::unit(c(0,-5,0,0), "cm")) +
  geom_confidence_tern(contour = TRUE,h = NULL, na.rm = FALSE, breaks = c(0.2, 0.5, 0.75),show.legend = NA, 
                       color='orange',size=0.5) +
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
             aes(pare5in, pare5out, pare5LR)) + 
  geom_point(color="gray10",size=1,shape=21) +
  theme_rgbw(base_size = 15) + 
  theme(plot.margin=grid::unit(c(0,0,0,-4), "cm")) +
  geom_confidence_tern(contour = TRUE,h = NULL, na.rm = FALSE, breaks = c(0.2, 0.5, 0.75),show.legend = NA, 
                       color='orange',size=0.5) +
  #theme(legend.justification=c(0,1), legend.position=c(0,1)) +
  #theme_gridsontop() + 
  #geom_confidence_tern(contour = TRUE, h = NULL, na.rm = FALSE, breaks = c(0.5, 0.9, 0.95), 
  #                    color='orange',linetype = "dashed") +
  guides(fill = guide_colorbar(order=1),color="none") + 
  labs(title = "Parental",
       x = "", xarrow = "Indegree", y = "", yarrow = "Outdegree",
       z = "", zarrow = "Local retension", size = 14) +
  theme(plot.title = element_text(hjust = 0.5, vjust= -8, size = 15)) + theme(legend.position = "none")
#labs(title= "Ternary Plot and Filled Contour",fill = "Value, V")

grid.arrange(b, b1, b2, b3, ncol = 2, nrow = 2)

#b1 + theme(plot.title = element_text(hjust = 0.5, size = 14)) + theme(legend.position = "none")

