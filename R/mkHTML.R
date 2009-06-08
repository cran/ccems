`mkHTML` <-
		function(smodels,showConstr=FALSE) {
	cat("\n ... making HTML file ... \n")
	.HTML.file=file(paste("results/", smodels[[1]]$msid,".htm",sep=""),"wt")
	cat("<html><h1>",smodels[[1]]$id,"Model Space</h1>",file = .HTML.file)
	cat("<h2>CPU time using ",smodels[[1]]$cpus," cpu(s) is ",file = .HTML.file)
	if (smodels[[1]]$totTime< 60) cat(format(smodels[[1]]$totTime,digits=3)," minutes",file = .HTML.file,sep="") else
	if (smodels[[1]]$totTime< 60*24) cat(format(smodels[[1]]$totTime/60,digits=3)," hours",file = .HTML.file,sep="") else 
		cat(format(smodels[[1]]$totTime/(60*24),digits=3)," days",file = .HTML.file,sep="")
	if (!is.null(smodels[[1]]$hists)) {
		cat("<h2>Model Space Size is ",smodels[[1]]$nMS," </h2><h2>|MS0| ... |MS5|: ",file = .HTML.file,sep="")
		upper=6
		for (i in 1:upper) cat(sprintf("%6d",smodels[[1]]$hists[i]),ifelse(i!=upper,", ","</h2>"),file = .HTML.file,sep="")
	}
	if (!is.null(smodels[[1]]$bestAics)) {
		cat("<h2>Best AICs: ",file = .HTML.file,sep="")
		upper=6
		for (i in 1:upper) cat(sprintf("%6s",format(smodels[[1]]$bestAics[i],digits=4)),ifelse(i!=upper,", ","</h2>"),file = .HTML.file,sep="")
	}
	cat("<p>",file = .HTML.file,sep="")
	cat("<p>",file = .HTML.file,sep="")
	if (!is.null(smodels[[1]]$ptype))
		cat("<h2>Parallelization method:  ",smodels[[1]]$ptype,"</h2>",file = .HTML.file,sep="")
	cat("<h2>Model space search method:  ",c("brute force","smart")[smodels[[1]]$smart+1],"</h2>",file = .HTML.file,sep="")
	cat("<h2>Model type:  ",ifelse(smodels[[1]]$TCC,"total concentration constraints","rational polynomial"),"</h2>",file = .HTML.file,sep="")
	repS=c("<TABLE BORDER=2 CELLPADDING=4>","<TR><TH>Model</TH><TH>Parameter</TH><TH>Initial Value</TH><TH>Optimal
					Value</TH><TH>Confidence Interval</TH></TR>")
	nMS=length(smodels)
# cat("\n nMS is ",nMS,"\n")
	for (i in 1:nMS) {
		g=smodels[[i]]
#		if (is.null(g$indx)) g$indx=0  # to avoid error if mkHTML is ever used on manually made and fitted model lists (bad idea since parallel at that point)
		nRows=dim(g$report)[1]
		biRcts=sum(g$W[row.names(g$report)[1],])-1 # number of binary reactions 
		# print("just before nOptParams in makeHTML")
		if ((g$nOptParams!=0) & (row.names(g$report)[1]!="p")  & 
          (row.names(g$report)[1]!="m1") & (length(grep("_",row.names(g$report)[1]))==0)
                      &  (length(grep("k",row.names(g$report)[1]))==0))
			repS=c(repS,sprintf("<TR><TH>%d %s.%d</TH><TD>%s</TD><TD>%4.3f^%d</TD><TD>%4.3f^%d</TD><TD>%s</TD></TR>",i,g$mid,g$indx,
			 row.names(g$report)[1],g$report[1,"initial"],biRcts,g$report[1,"final"]^(1/biRcts),biRcts,g$report[1,"confidenceInterval"])) else
			repS=c(repS,sprintf("<TR><TH>%d %s.%d</TH><TD>%s</TD><TD>%4.3f</TD><TD>%4.3f</TD><TD>%s</TD></TR>",i,g$mid,g$indx,
						row.names(g$report)[1],g$report[1,"initial"], g$report[1,"final"],g$report[1,"confidenceInterval"])) 
		for (j in 2:nRows)		{
			biRcts=sum(g$W[row.names(g$report)[j],])-1 # number of binary reactions 
			if (j <= ifelse(showConstr,g$nZ+1,g$nOptParams))  
				if ((row.names(g$report)[j]!="p")  & (row.names(g$report)[j]!="m1")  & (length(grep("_",row.names(g$report)[j]))==0) 
                         & (length(grep("k",row.names(g$report)[j]))==0))
					repS=c(repS,sprintf("<TR><TD>&nbsp;</TD><TD>%s</TD><TD>%4.3f^%d</TD><TD>%4.3f^%d</TD><TD>%s</TD></TR>",
				row.names(g$report)[j],g$report[j,"initial"],biRcts,g$report[j,"final"]^(1/biRcts),biRcts,g$report[j,"confidenceInterval"])) else
					repS=c(repS,sprintf("<TR><TD>&nbsp;</TD><TD>%s</TD><TD>%4.3f</TD><TD>%4.3f</TD><TD>%s</TD></TR>",
									row.names(g$report)[j],g$report[j,"initial"], g$report[j,"final"],g$report[j,"confidenceInterval"]))
			if (j > ifelse(g$activity,2*g$nZ,g$nZ)+2) 
				repS=c(repS,sprintf("<TR><TD>&nbsp;</TD><TD>%s</TD><TD>%4.3f</TD><TD>%4.3f</TD><TD>%s</TD></TR>",
								row.names(g$report)[j],g$report[j,"initial"], g$report[j,"final"],g$report[j,"confidenceInterval"]))
		} 
	} # loop on i
	repS=c(repS,"</TABLE>")
	writeLines(repS,con=.HTML.file)
	cat("\n</html>", append = TRUE,file = .HTML.file)
	close(.HTML.file)
} # makeHTML

