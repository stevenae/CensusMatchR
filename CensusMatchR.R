######
# Title: CensusMatchR
# Version: 0.1
######
# Author:	Steven Ellis
#			elliss@google.com
######
# LICENSE (GNU GENERAL PUBLIC LICENSE 2)
# Copyright 2013 Google Inc.
######
# AUTHORS' CAUTION:
#   This code is for RESEARCH PURPOSES only and has no warranty. 
#   It almost certainly contains large and small errors. 
#   Evaluate it thoroughly for your own purposes. Read the code!
######
# REFERENCES:
#   Identifying Surrogate Geographic Research Regions with Advanced Exact Test Statistics (2013).
#   AMA ART Forum, Chicago, IL.
######
# DATA:
#   http://www2.census.gov/acs2010_1yr/pums/csv_hus.zip
#   http://wonder.cdc.gov/wonder/sci_data/datasets/zipctyA.zip
#   http://wonder.cdc.gov/wonder/sci_data/datasets/zipctyB.zip
#   http://www.mcc.co.mercer.pa.us/dps/state_fips_code_listing.htm [saved as three-column CSV, loaded as StateFips]
######


# Example code, matching Chicago, IL to points around the nation.
# Matches on Household Income, Household Size, Household Ownership, and Units in Structure.
if(FALSE) {
	ZipDir <- "/small_zip_files/" # place both all_zips.csv from zipctyA and zipctyB here, as they arrive or split these into smaller files
	StatesDir <- "~/by_state_hh/"
	StateFipsFile <- "~/state_fips.csv"
	PumaCodeMapFile <- "~/new_puma.csv"
	RegionFilenames <- c("ss10hil.csv")	# you may specify multiple states
	RegionPumas <- seq(03501, 03519, 1)
	MatchVars <- c("HINCP","NOC","TEN","BLD")
	Matches <- RunModel(1000) # use lower sample size to test, higher for greater accuracy
	ZipCodes <- MatchesToZips(Matches)
}

# Abstracted function to process surrogate regions and return summary statistics.
RunModel <- function(PrepLength) {
	require(MASS)
	require(crossmatch)
	require(foreach)
	require(utils)
	require(doMC)

	RegionFilenames <- paste(StatesDir, RegionFilenames, sep='')
	RegionList <- lapply(RegionFilenames, function(x) {
		DFIntoPumaList(read.csv(x, strip.white = TRUE))
	})
	RegionList <- Reduce(c,RegionList)
	RegionList <- RegionList[which(names(RegionList) %in% RegionPumas)]

	Timing <- list()
	CompletedMatches <- CompleteMatches(PrepLength, RegionList)
	return(CompletedMatches)
}

# Takes output of RunModel() and Matches PUMA codes to Zip codes for real-world use and mapping.
MatchesToZips <- function(CompletedMatches) {
	SignificantMatchByRegionPuma <- lapply(CompletedMatches, function(RegionPuma){
		lapply(RegionPuma, function(State) {
			State <- State[!is.null(State)]
			Matches <- lapply(State, function(x) {
				if(is.null(x)){
					return(NA)
				}
				else if (x$approxpval < .05) {
					x$a1
				}
				else {
					return(NA)
				}
			})
			Matches[!is.na(Matches)]
		})
	})

	UnZipped <- Reduce(c,Reduce(c, SignificantMatchByRegionPuma))
	
	StateFips <- read.csv(StateFipsFile, strip.white = TRUE)
	PumaCodeMap <- read.csv(PumaCodeMapFile, strip.white = TRUE)

	Zipped <- c() 
	for (x in list.files(ZipDir)) {
		print(paste("loading file ",x,"...", sep=""))
		ZipFile <- read.csv(paste(ZipDir, x,sep='/'), strip.white = TRUE)
		Zipped <- c(Zipped,ZipFileZipCodes(ZipFile, UnZipped, StateFips, PumaCodeMap))
	}
	Zipped <- as.list(rev(sort(unlist(Zipped))))
	return(Zipped)
}

# Worker function for MatchesToZips(), to allow for flexible processing of individual zip code files.
ZipFileZipCodes <- function(ZipFile, PumaList, StateFips, PumaCodeMap) {
	print(paste("loaded, analyzing",as.character(nrow(ZipFile)),"records"))
	pb <- txtProgressBar(min = 0, max = length(PumaList), style = 3)
	registerDoMC()

	SignificantZips <- foreach (i = icount(length(PumaList)), .combine='c') %dopar% {
		setTxtProgressBar(pb, i)
		CurrCode <- PumaList[i]
		StFip <- PumaCodeMap$B[which(PumaCodeMap$D == names(CurrCode))]
		StateText <- as.vector(na.omit(StateFips$State.Abbreviation[which(StateFips$FIPS.Code == StFip)]))
		CoFip <- na.omit(unique(PumaCodeMap$E[which(PumaCodeMap$D == names(CurrCode))]))

		StateRows <- which(ZipFile[,7] == StateText)
		FipsRows <- which(ZipFile[,8] == CoFip)
		
		NewZip <- unique(ZipFile[intersect(StateRows, FipsRows),1])
		if (length(NewZip)) {
			for (i in seq(1,length(NewZip))) {
				if(nchar(NewZip[i]) == 4) {
					NewZip[i] <- paste("0",NewZip[i],sep="")
				}
			}
			CurrCode <- rep(CurrCode,length(NewZip))
			names(CurrCode) <- NewZip
			CurrCode
		}
		else {
			NULL
		}
	}
	close(pb)
	return(SignificantZips)
}

# Retrieving test statistics for single potential surrogate region.
CompleteMatches<- function(PrepLength, RegionList) {
    CompletedMatches <- NULL
    CompletedMatches <- lapply(RegionList, function(RegionPuma) {
		print("new Puma running...")
        lapply(list.files(StatesDir, pattern = 'csv'), function(x) {
            tryCatch(ComparePumas <- DFIntoPumaList(read.csv(paste(StatesDir, x, sep = ''), strip.white = TRUE)),silent = TRUE,
                error = function(e) {
                NULL
                }, finally = NA)
            NextPuma <- lapply(ComparePumas, function(ActPuma) {
                print(PrepLength)
                matched <- NULL
                tryCatch(Timing[[as.character(x)]] <- system.time(prepped <- PrepareForCrossmatch(FixList(RegionPuma), FixList(ActPuma), PrepLength)),silent = TRUE,
                    error = function(e) {
                    NULL
                    }, finally = NA)
                tryCatch(print(as.character(x)),silent = TRUE,
                    error = function(e) {
                    NULL
                    }, finally = NA)
                tryCatch(print(Timing[[as.character(x)]]),silent = TRUE,
                    error = function(e) {
                    NULL
                    }, finally = NA)
                tryCatch(print(system.time(matched <- crossmatchtest(prepped[[1]], prepped[[2]]))),silent = TRUE,
                    error = function(e) {
                    NULL
                    }, finally = NA)
                return(matched)
            })
        })
    })
    return(CompletedMatches)
}

# Preparation of covariance matrices for matching.
PrepareForCrossmatch <- function(RegionList, CompareList, PrepLength ) {
	CompareListSub <- lapply(CompareList, function(xx) {
        sample(as.vector(unlist(xx)), PrepLength, replace = T)
    })
	TestCompare<-do.call(rbind, CompareListSub)
	RegionCompareListSub <- lapply(RegionList, function(xx) {
        sample(as.vector(unlist(xx)), PrepLength, replace = T)
    })
	RegionCompare<-do.call(rbind, RegionCompareListSub)
	z <- c(rep(0,PrepLength), rep(1,PrepLength))
	X=t(cbind(TestCompare, RegionCompare))
	X<-as.matrix(X)
	n<-dim(X)[1]
	k<-dim(X)[2]
	for (j in 1:k) X[,j]<-rank(X[,j])
	cv<-cov(X)
	vuntied<-var(1:n)
	rat<-sqrt(vuntied/diag(cv))
	cv<-diag(rat)%*%cv%*%diag(rat)
	out<-matrix(NA,n,n)
	icov<-ginv(cv)
	for (i in 1:n) {
		out[i,]<-mahalanobis(X,X[i,],icov,inverted=TRUE)
	}
	dis<-out
	return(list(z, dis))
}

# Data munging.
FixList <- function(RawList) {
	RawList <- lapply(RawList, function(RL) {
		as.vector(na.exclude(RL[which(RL>0)]))
	})

	RawLengths <- lapply(RawList, function(RL) {
		length(RL)
	})
	MaxLength <- max(unlist(RawLengths))
	RawList <- lapply(RawList, function(RL) {
		if (length(RL) < MaxLength) {
			c(RL, sample(RL, abs(length(RL) - MaxLength), replace = TRUE))
		}
		else {
			RL
		}
	})
	return(RawList)
}

# Select only variables of interest from a given PUMA's data.
DFIntoPumaList <- function(df) {
	InterestList <- lapply(unique(df$PUMA), function(PumaCode) {
		tryCatch(GetDFInfo(df[which(df$PUMA == PumaCode),]),silent = TRUE,
            error = function(e) {
            return(NULL)
            }, finally = NA)
	})
	names(InterestList) <- unique(df$PUMA)
	return(InterestList)
}

# Factored logic for DFIntoPumaList().
GetDFInfo <- function(df) {
	ReturnList<-lapply(MatchVars, function(MatchVar) {
		df[,MatchVar]
	})
	names(ReturnList) <- MatchVars
	return(ReturnList)
}
