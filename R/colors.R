library(RColorBrewer)

colors=c("#E04C5C", "#FB894B", "#E7DA36", "#7DAF4C", "#187A51",
         "#5EA4A2", "#23AECE", "#3D3C4E", "#4D1836", "#C51B7D",
         "#E9A3C9",  "#B35806", "#F1A340", "#FEE08B", "#D9EF8B",
         "#91CF60", "#C7EAE5", "#5AB4AC", "#01665E", "#E7D4E8",
         "#AF8DC3", "#762A83")

# arm.color <- c('maternal Tdap - infant ap' = "#E04C5C", 'maternal Tdap - infant wp' = "#7DAF4C", 'No maternal vaccination - infant wp' = "#23AECE")
# arm.short.color <- c('matplus_iap' = "#E04C5C", 'matplus_iwp' = "#7DAF4C", 'matminus_iwp' = "#23AECE")
arm.reserve.color <- c("#E04C5C", "#7DAF4C", "#23AECE", "#FB894B", "#E7DA36", "#187A51",
                       "#5EA4A2", "#3D3C4E", "#4D1836", "#C51B7D", "#E9A3C9", "#B35806", "#F1A340", "#FEE08B", "#D9EF8B",
                       "#91CF60", "#C7EAE5", "#5AB4AC", "#01665E", "#E7D4E8", "#AF8DC3", "#762A83")

# arm.boxplot.color <- c("#E04C5C", "#86B3B6", "#23AECE", "#E9A3C9")
# arm.crossbar <- c("#E04C5C","#E04C5C","#86B3B6","#86B3B6", "#23AECE", "#23AECE")
#
# arm.apwp.color <- c('maternal Tdap - infant ap' = "#E04C5C", 'maternal Tdap - infant wp' = "#7DAF4C")
# arm.apwps.color <- c('mat+_iap' = "#E04C5C",'mat+_iwp' = "#7DAF4C")
# arm.apwps.color <- c('mat+_iap' = "dodgerblue2",'mat+_iwp' = "forestgreen")

arm.apwps.color <- c('aP' = "dodgerblue2",'wP' = "forestgreen")

arm.matpm.color <- c('mTT' = "salmon",'mTdaP' = "skyblue")

# arm.matpm.color <- c('maternal Tdap - infant wp' = "#7DAF4C", 'No maternal vaccination - infant wp' = "#23AECE")
# arm.matpms.color <- c('Tdap+' = "#7DAF4C", 'Tdap-' = "#23AECE")
# arm.matpm.boxplot.color <- c("#7DAF4C", "#23AECE")
# arm.matpm.boxplot.crossbar <- c("#7DAF4C", "#7DAF4C", "#23AECE", "#23AECE")

sex.color <- c('Female' = "#F1A340", 'Male' = "#91CF60")

# visit.color <- c('vaccinated' = "white", 'boosted' = "lightgrey")
visit.color.arm <- c('prevaccinated' = "#5AB4AC", 'vaccinated' = "#C51B7D")
# visit.shape <- c('vaccinated' = 21, 'boosted' = 16)
# visit.boxplot.color <- c("#4D1836", "#91CF60")
# visit.crossbar <- c("#4D1836", "#4D1836", "#91CF60", "#91CF60")

clust.colors <- c("#E69F00", "#56B4E9", "#009E73", "#D31E00", "#0072B2", "#D55E00", "#CC79A7")
tangle.clust.colors <- c('1'="#E69F00", '2'="#56B4E9", '3'="#009E73", '4'="#D31E00", '5'="#0072B2", '6'="#D55E00", '7'="#CC79A7")
clust.color.umap <- c('1' = "#E69F00", '2' = "#56B4E9", '3' = "#009E73", '4' = "#D31E00", '5' = "#0072B2", '6' = "#D55E00",'7' = "#CC79A7")

clust.colors.8 <- c("#E69F00", "#56B4E9", "#009E73", "#D31E00", "#762A83", "#0072B2", "#D55E00", "#CC79A7")
tangle.clust.colors.8 <- c('1'="#E69F00", '2'="#56B4E9", '3'="#009E73", '4'="#D31E00", '5'="#762A83", '6'="#0072B2", '7'="#D55E00", '8'="#CC79A7")
clust.color.umap.8 <- c('1' = "#E69F00", '2' = "#56B4E9", '3' = "#009E73", '4' = "#D31E00", '5'="#762A83", '6' = "#0072B2", '7' = "#D55E00",'8' = "#CC79A7")


# tangle.clust.colors <- c('1'="red",'2'="darkolivegreen",'3'="deepskyblue",'4'="deeppink2",'5'="darkorchid3",'6'="darkorange2",'7'="cyan4")

anno_colors <- list(Arm = arm.apwps.color)

# arm.s.color <- c("#7DAF4C", "#23AECE")
# names(arm.s.color) <- c("Tdap+", "Tdap-")

#Antigen coloring for dendrogram labels
dendrogram.antigen.colors <- c("PT"="orangered","TT"="blue3","PRN"="orange2","DT"="green3","FHA"="gold")


#Feature coloring for dendrogram labels
dendrogram.feature.colors <- c("IgG"="darkgrey","IgG1"="orangered","IgG2"="lightblue3","IgG3"="orange2","IgG4"="blue3","FcgR2a"="green3","FcgR3b"="green3","ADCP"="purple","ADNP"="purple","ADCD"="pink")

##arm.color <- c('Non-Pregnant 15 mcg' = "#E9A3C9", 'Pregnant 15 mg' = "#7DAF4C", 'Pregnant 30 mcg' = "#187A51")

##arm.color <- c('Non-Pregnant 15 mcg' = "#D37295", 'Pregnant 15 mg' = "#86B3B6", 'Pregnant 30 mcg' = "#499894")

