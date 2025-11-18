# manhattan plot per treat, plus the cmh clear model after
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(viridis)) install.packages('viridis')
if (!require(ggpubr)) install.packages('ggpubr')
if (!require(scales)) install.packages('scales')
if (!require(ggnewscale)) install.packages('ggnewscale')
if (!require(svglite)) install.packages('svglite')
if (!require(extrafont)) install.packages('extrafont')

# create plots for THB, LNHB, TNHH, LNHH in that order
treat <- c("THB",  "LHB", "TNHB", "LNHB", "THH", "LHH", "TNHH", "LNHH" )

# read in model plots -------------------------------------------------------------------------------------
# read in data
setwd("~/mol_eco_2024/2.5.3_sel_scans")
mod_cmh <- fread(sep = "\t", paste0("model_both_clear_cmh_num_snps_soil_sep.tsv"))

# set cols
#"#0071c1ff", "#538234ff") # colors from pub 2024 nat com
soil_colors <- c(Limestone = "#0071c1ff", Tuff = "#538234ff")
herb_shapes <- c(H = 24, NH = 21)

# set soils so plots nicer and numbered groups to letters
mod_cmh$soil <- ifelse(mod_cmh$soil == "L", "Limestone", "Tuff")
mod_cmh$bee <- ifelse(mod_cmh$bee == "H", "Hand", "Bee")
mod_cmh$group_let <-letters[mod_cmh$soil_group]

mod_cmh$let <- ifelse(mod_cmh$soil == "Limestone", toupper(mod_cmh$group_let), mod_cmh$group_let)

# one that giacomo, emi and corrine suggested
cmh_t <- mod_cmh[soil == "Tuff"]
cmh_l <- mod_cmh[soil == "Limestone"]

n_max <- 1.15*max(mod_cmh$mean)



# make subsets for easier annotation
#t_sub <- cmh_t[herb == "NH" & bee == "Hand"]
#l_sub <- cmh_l[herb == "NH" & bee == "Bee"]


prop_snps <- ggarrange(
	ggplot(cmh_t, aes(x = bee, fill = soil, shape = herb)) +
	    geom_errorbar(aes(y = mean, ymin= lcl, ymax= ucl), width=.5, position= position_dodge(width= 1)) +
	    geom_point(aes(y = mean, shape = herb), size= 7, position= position_dodge(width= 1)) +
	    geom_text(aes(y = ucl*0.96, label = let, family= "Arial"), color = "black", size = 5, vjust = -1.5, position= position_dodge(width= 1)) +
	    scale_fill_manual(values=soil_colors) +
	    scale_shape_manual(values=herb_shapes) +
	    scale_y_continuous(limits = c(0, n_max)) + 
	    labs(x= NULL, y= NULL) +
	    theme_bw(base_size = 22) +
	    theme(text = element_text(family = "Arial")),
	ggplot(cmh_l, aes(x = bee, fill = soil, shape = herb)) +
	    geom_errorbar(aes(y = mean, ymin= lcl, ymax= ucl), width=.5, position= position_dodge(width= 1)) +
	    geom_point(aes(y = mean, shape = herb), size= 7, position= position_dodge(width= 1)) +
	    geom_text(aes(y = ucl*0.96, label = let, family= "Arial"), color = "black", size = 5, vjust = -1.5, position= position_dodge(width= 1)) +
	    scale_fill_manual(values=soil_colors) +
	    scale_shape_manual(values=herb_shapes) +
	    labs(x= NULL, y= NULL) +
			scale_y_continuous(limits = c(0, n_max)) + 	
	    theme_bw(base_size = 22) +
	    theme(text = element_text(family = "Arial")),
ncol= 2, nrow = 1, labels = c("A","B"),  font.label = (family = "Arial"), common.legend = TRUE, legend = "right") # leave legend so easier to change in inkscape

prop_snps <- annotate_figure(prop_snps, left = text_grob("Proportion of Evolved SNPs", 
		color = "black", size = 14, rot = 90, family = "Arial"))
                    #bottom = text_grob("Pollination Treatment", color = "black", face = "bold", size = 14))

# to read in cmh --------------------------------------------------------------------------------------------------------------------------
setwd("~/mol_eco_2024/2.5.3_sel_scans")
clear_total <- fread(sep = "\t", paste0("clear_and_cmh_all_treat_pruned.tsv"))

# look up table
x1 <- c(1,2,3,4,5,6,7,8,9,10)
x2 <- c("A01","A02","A03","A04","A05","A06","A07","A08","A09","A10")
look <- data.table(x1,x2)
setkey(look, x2)

# merge to get numerical chrome
setkey(clear_total, chr)
clear_total <- clear_total[look]
setnames(clear_total, c("x1"), c("num_chr"), skip_absent=TRUE)	

# create ID for easier checking if there
clear_total[, id := paste(num_chr, pos, sep = "_")]


# for plotting later
odd <- seq(from = 1, to = 10, by = 2)
even <- seq(from = 2, to = 10, by = 2)

# get max for axis
y_max <- max(-log10(clear_total[treatment %in% treat]$fdr))

for(i in 1:(length(treat))) {
	subset = treat[i]
	print(paste0("loading files for subset ", subset, " - - - - - - - - - - - - - - - - - - - - - - - - -"))

	# read in total cmh with updated ne
	cmh <- clear_total[treatment == subset] 
	setkey(cmh, num_chr)

	# add get max pos per chr
  	max_pos <- cmh[,  max(pos), by = num_chr]
  	setnames(max_pos, "V1", "max_bp")

  	# get max from prev chrom to add
  	max_pos[, bp_add := lag(cumsum(max_bp), default = 0)] 
  	setkey(max_pos, num_chr)

  	# merge to cmh and add the max from the prev chrom 
  	cmh <- cmh[max_pos]
  	cmh[, bp_cum := pos + bp_add]

	# get mean axis to center on
	axis_set <- cmh[, mean(bp_cum), by=c("num_chr")] # mean or median?
	setkey(axis_set, "num_chr")

	# create rectangles to shade
  	rect_min <- cmh[, min(bp_cum), by=c("num_chr")]
	rect_max <- cmh[, max(bp_cum), by=c("num_chr")]
	setkey(rect_min, "num_chr")
	setkey(rect_max, "num_chr")
	rect <- rect_min[rect_max]
	rect$num_chr <- as.factor(rect$num_chr)

	# set color, dark grey if odd, light grey if even
	cmh$col <- ifelse(cmh$num_chr %in% odd, 
                       "A", 
                       ifelse(cmh$num_chr %in% even, "B", NA))

	# yellow if in gwas, purple if just under sel
	cmh$col <- ifelse(cmh$sig == TRUE, "D", cmh$col)

	# load magma cols
	m <- magma(3)
	group.colors <- c(A = "grey70", B = "grey90", C =m[2], D = m[1])

  	###########################################################################################################
  	# plot 
	p1 <- ggplot() +
		geom_point(data = cmh, aes(x = bp_cum, y = -log10(fdr), color = col, ), alpha = 0.8) +
		scale_color_manual(values=group.colors) +

		# create rectangles to shade over, can hide if we want, can also play with alpha
		new_scale_colour() +
		geom_rect(data = rect, aes(xmin=rect$V1, xmax=rect$i.V1, ymin = 0.6, ymax = Inf, fill = rect$num_chr), alpha = .04) + 
		scale_fill_manual(values = c(rep_len(c("white", "grey10"), length(unique(rect$num_chr))))) +

		# create axises, may want to hide legend title as it is repetitive
		scale_x_continuous(label = axis_set$num_chr, breaks = axis_set$V1, , expand = c(0.01, 0)) +
		scale_y_continuous(limits = c(0, y_max*1.05), expand = c(0.01, 0.01)) +
		#scale_y_continuous(expand = c(0,0), limits = c(1, max(res_1$pleio)*1.2)) +
		labs(x = NULL, y = NULL) + 
		theme_minimal(base_size = 16) +
		theme(legend.position = "none",
			panel.grid = element_blank(),
			axis.line = element_line(color="black", linewidth = 0.25),
        axis.ticks = element_line(size = 0.25, color="black"),
        text = element_text(family = "TT Arial"))

	# save for adding together
	assign(paste0("p1_", subset), p1) # save p1 for plotting of total later
	print(paste0("Saved first plot for subset ", subset))	

}

# loop over saved plots to add together ------------------------------------------------------------------------------------------------
rm(p1)
p_list <- paste("p1", treat, sep = "_")

# automating the body of code below
all_plots <- noquote(paste(p_list, collapse = ", "))

# add the letters in front of the
l <- LETTERS[seq( from = 3, to = 2+length(treat) )]
#t_lab <- paste(l, treat, sep = "   ")
#plot_labs <- c(paste0(sprintf("'%s'", t_lab), collapse = ", "))

let_labs <- c(paste0(sprintf("'%s'", l), collapse = ","))

# arrange and paste with labels
str2 <-paste0("totaL_pleio_plot <- ggarrange(", all_plots, ", labels = c(",let_labs,"), font.label = (family = 'Arial'), hjust = -0.9,
	ncol = 2, nrow = length(treat)/2, common.legend = TRUE, legend='none')")
eval(parse(text = str2))

# add axis labels
totaL_pleio_plot <- annotate_figure(totaL_pleio_plot, 
			left = text_grob(expression(paste("-log"[10], "(FDR(", italic("p"), "))")), 
				color = "black", size = 14, rot = 90, family = "Arial"),
            bottom = text_grob(paste0("Chromosomes"), color = "black", size = 14, family = "Arial"))

# arrange the model and the manhattan plots 
model_man_plot <- ggarrange(
	prop_snps,
	totaL_pleio_plot,
	ncol= 1, nrow = 2, labels = NULL, common.legend = TRUE, legend="none", 
	widths = c(1, 1), heights = c(1, 1.33333333)) # adjust the relative width and height of plots, I think it only works for cols tho
model_man_plot


resolution <- 144
svglite(paste0("both_clear_cmh_model_manhattan_all_treat.svg", sep=""), width = 1440/resolution, height = 1280/resolution)

model_man_plot

dev.off()

# save as png
resolution <- 144
ggsave(filename = "both_clear_cmh_model_manhattan_all_treat.png", plot = model_man_plot, 
	width = 1440/resolution, height = 1280/resolution, units = "in", dpi= "retina", bg = "white")