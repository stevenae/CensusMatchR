library(ggplot2)
library(zipcode)
library(RColorBrewer)
data(zipcode)

newzips <- mclapply(names(Zipped), function(x) {
	if(nchar(x) == 4) {
		x <- paste(x,"0",sep="")
	}
	})

zipcode$plotter <- 0
zipcode$plotter[which(zipcode$zip %in% names(Zipped))] <- 1
zipcode$plotter <- as.factor(zipcode$plotter)
levels(zipcode$plotter) <- c("Insignificant","Matched")

pdf(width=20, height=10)
g = ggplot(data=zipcode) + geom_point(aes(x=longitude, y=latitude, colour=rev(plotter), position="jitter"))
g = g + scale_colour_manual(values = brewer.pal(8,"Blues")[c(2,8)], name="Algorithm Result")
# simplify display and limit to the "lower 48"
g = g + theme_bw() + scale_x_continuous(limits = c(-125,-66), breaks = NA)
g = g + scale_y_continuous(limits = c(25,50), breaks = NA)
# don't need axis labels
g = g + labs(x=NULL, y=NULL)
ggsave(filename="~/ama_art_map.pdf", g)
dev.off()