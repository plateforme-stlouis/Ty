

library(readr)
library(ggplot2)
library(doBy)

dat_full <- read_delim("inFile.txt",
                       "\t", escape_double = FALSE, col_names = FALSE,
                        trim_ws = TRUE)

colnames(dat_full) <- c("qname"
                      , "flag"
                      , "rname"
                      , "pos"
                      , "mapq"
                      , "cigar"
                      , "rnext"
                      , "pnext"
                      , "tlen"
                      , "seq"
                      , "qual"
                        )

is.GoodQuality <- Reduce('&&', dat_full$mapq > 0)
if (! is.GoodQuality) {
    distribution.quality <- data.frame(chr = dat_full$rname
                                     , mapq = dat_full$mapq)
    distribution.quality$num <- 1
    distribution.quality <- summaryBy(num ~ mapq + chr
                                    , data = distribution.quality
                                    , FUN = length
                                    , keep.names = TRUE)
    distribution.quality.all <- summaryBy(num ~ mapq
                                        , data = distribution.quality
                                        , FUN = sum
                                        , keep.names = TRUE)
}

flags <- c(83, 99)
dat_filter1 <- subset(dat_full, flag %in% flags)
dat_filter2 <- data.frame(flag   =dat_filter1$flag
                        , chr    =dat_filter1$rname
                        , posBeg =dat_filter1$pos
                        , posEnd =dat_filter1$pnext
                        , posDiff   =rep(0, nrow(dat_filter1))
                        , pos   =rep(0, nrow(dat_filter1))
                        , tlen   =dat_filter1$tlen
                        , len    =rep(0, nrow(dat_filter1))
                        #, orient =rep("+", nrow(dat_filter1))
                        , orient =rep(+1, nrow(dat_filter1))
                        , dup    =rep(0, nrow(dat_filter1))
                        , weight    =rep(0, nrow(dat_filter1))
                        , seq    =dat_filter1$seq
                        , stringsAsFactors=FALSE
                          )
dat_filter2[dat_filter2$flag==83, "orient"] <- -1
dat_filter2[, "len"] <- apply(matrix(dat_filter2$seq), 1, function(x){nchar(as.character(x))})
dat_filter2[, "posDiff"] <- apply(matrix(c(dat_filter2$posBeg, dat_filter2$posEnd), ncol=2)
                                    , 1
                                    , function(x){x[2] - x[1]})
dat_filter2[, "pos"] <- apply(matrix(c(dat_filter2$orient, dat_filter2$posBeg, dat_filter2$posEnd)
                                   , ncol=3)
                            , 1
                            , function(x){
                                if (x[1]>0) {
                                    x[2]
                                } else {
                                    x[3]
                                }
                            })
dat <- dat_filter2

udat <- unique(dat_filter2)
ddat <- dat[duplicated(dat), ]

isAll_same_len <- Reduce('&&', dat[, "len"] == dat[1, "len"])
isAll_same_orient <- Reduce('&&', sign(dat[, "tlen"]) == sign(dat[, "orient"]))

for (i in 1:nrow(ddat)) {
    ndup <- 0
    ind <- NULL
    for (j in 1:nrow(dat)) {
        if (Reduce('&&', ddat[i, ] == dat[j, ])) {
            ind <- c(ind, j)
            ndup <- ndup + 1
        }
    }
    dat[ind, "dup"] <- ndup
}
## isAll_own_duplicate <- Reduce('&&', dat[, "dup"] == 1)
isAll_own_duplicate <- nrow(ddat) == 0

dat1 <- summaryBy(dup ~ chr + posBeg + orient + seq
                , FUN=length
                , dat
                , keep.names=TRUE)
dat2 <- summaryBy(dup ~ chr + posBeg + orient
                , FUN=length
                , dat
                , keep.names=TRUE)

res <- summaryBy(weight ~ chr + pos + orient
               , FUN=length
               , udat
               , keep.names=TRUE)
