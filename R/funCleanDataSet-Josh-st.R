

#'
#' Read one file and return the data.frame
#'
#' @param mypath the path of the file
#' @return the data.frame with the columns: chr pos orient weight
#'
#' Need to be a bit more cleanned !! :-)
#' (some stuff are not useful anymore)
#'
#' @examples
#'
#' AB15 <- system.file("extdata", "AB15.light.txt", package = "Ty")
#' AB16 <- system.file("extdata", "AB16.light.txt", package = "Ty")
#'
#' r15 <- GetWeight(AB15)
#' r16 <- GetWeight(AB16)
#'
#'
#' @import readr doBy
#'
#' @export
GetWeight <- function(mypath) {
    dat_full <- read_delim(mypath
                         , " "
                         , escape_double = FALSE
                         , col_names     = FALSE
                         , trim_ws       = TRUE)

    print(paste("Read done. From", mypath))

    colnames(dat_full) <- c("qname"
                          , "flag"
                          , "rname"
                          , "pos"
                            #, "mapq"
                            #, "cigar"
                            #, "rnext"
                          , "pnext"
                            #, "tlen"
                          , "seq"
                            #, "qual"
                            )


    print("Cleaning up...")

    flags <- c(83, 99)
    dat_filter1 <- subset(dat_full, flag %in% flags)
    dat_filter2 <- data.frame(flag             =dat_filter1$flag
                            , chr              =dat_filter1$rname
                            , posBeg           =dat_filter1$pos
                            , posEnd           =dat_filter1$pnext
                            , posDiff          =rep(0, nrow(dat_filter1))
                            , pos              =rep(0, nrow(dat_filter1))
                              #, tlen          =dat_filter1$tlen
                            , len              =rep(0, nrow(dat_filter1))
                              #, orient        =rep("+", nrow(dat_filter1))
                            , orient           =rep(+1, nrow(dat_filter1))
                            , dup              =rep(0, nrow(dat_filter1))
                            , weight           =rep(0, nrow(dat_filter1))
                            , seq              =dat_filter1$seq
                            , stringsAsFactors =FALSE
                              )
    dat_filter2[dat_filter2$flag==83, "orient"] <- -1
    ## dat_filter2[, "len"] <- apply(matrix(dat_filter2$seq), 1, function(x){nchar(as.character(x))})
    ## dat_filter2[, "posDiff"] <- apply(matrix(c(dat_filter2$posBeg, dat_filter2$posEnd), ncol=2)
    ##                                 , 1
    ##                                 , function(x){x[2] - x[1]})
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

    print("Computing weight...")

    udat <- unique(dat_filter2)

    res <- summaryBy(weight ~ chr + pos + orient
                   , FUN=length
                   , udat
                   , keep.names=TRUE)

    print("Weight computed.")


    # Hum? TODO: fix all that !!
    res$chrX <- sapply(res$chr, Ref2Str)


    res$From <- rep("Case")

    # Ugly conversion
    res$pos <- as.integer(res$pos)
    res$weight <- as.integer(res$weight)
    res$chr <- as.factor(res$chr)
    res$chrX <- as.factor(res$chrX)
    res$From <- as.factor(res$From)

    res
}


#'
#' Merge only one chromosome from two data.frames
#'       (?not useful?)
#'
#' The data.frames come from the function GetWeight
#'
#' @param df1 first data.frame
#' @param df2 second data.frame
#' @param kro name of one chromosome
#'
#' @return the merged data.frame
#'
#'
#' df1 comes from AB15
#' df2 comes from AB16
#'
#' kro is "chrI" or "ref|NC_001133|"
#'
#' @export
MergeOnly <- function(df1, df2, kro = "chrXX") {
    s1 <- subset(df1, df1$chr == kro)
    s2 <- subset(df2, df2$chr == kro)

    r <- merge(s1, s2, by = c("pos", "chr"))

    r
}

#'
#' Bind by row only one chromosome from two data.frames
#'
#' The data.frames come from the function GetWeight
#'
#' @param df1 first data.frame
#' @param df2 second data.frame
#' @param kro name of one chromosome
#' @param name1 associates a name to df1
#' @param name2 associates a name to df2
#'
#' @return the binded data.frame by row
#'
#' Note:
#'
#' df1 comes from AB15
#' df2 comes from AB16
#'
#' kro is "chrI" or "ref|NC_001133|"
#'
#' @export
rBindOnly.Signed <- function(df1, df2, kro = "chrXX", name1 = "One", name2 = "Two") {
    s1 <- subset(df1, df1$chr == kro)
    s2 <- subset(df2, df2$chr == kro)

    s2$weight <- - s2$weight

    s1$name <- rep(name1, nrow(s1))
    s2$name <- rep(name2, nrow(s2))

    r <- rbind(s1, s2)

    r$pos.OrdFac <- factor(r$pos, ordered=TRUE)

    r
}

#'
#' Randomly select n rows from data.frame
#'
#' @param df the data.frame
#' @param n the number of randomly selected rows (default: 20)
#' @return a subset and tiny data.frame
#'
#' @export
GetTinySample <- function(df, n = 20) {
    if (n > nrow(df)) {
        n <- nrow(df)
    }
    samp <- df[sample(1 : nrow(df), n, replace = FALSE),]

    samp
}


#'
#' Table that binds weird Identification to pronounceable name
#'
#' then used by a useful? function
#'
#' note: stringAsFactors returns Error. Hum?
#' therefore, apply as.character to the right columns
#'
#' @export
chr.ref.NC <- read.table(header = TRUE, text = "
chr_num  Chromosome.X     Genbank_ID RefSeq_ID  Length   chrX
1       'Chromosome I'    BK006935.2 NC_001133  230218   chrI
2       'Chromosome II'   BK006936.2 NC_001134  813184   chrII
3       'Chromosome III'  BK006937.2 NC_001135  316620   chrIII
4       'Chromosome IV'   BK006938.2 NC_001136 1531933   chrIV
5       'Chromosome V'    BK006939.2 NC_001137  576874   chrV
6       'Chromosome VI'   BK006940.2 NC_001138  270161   chrVI
7       'Chromosome VII'  BK006941.2 NC_001139 1090940   chrVII
8       'Chromosome VIII' BK006934.2 NC_001140  562643   chrVIII
9       'Chromosome IX'   BK006942.2 NC_001141  439888   chrIX
10      'Chromosome X'    BK006943.2 NC_001142  745751   chrX
11      'Chromosome XI'   BK006944.2 NC_001143  666816   chrXI
12      'Chromosome XII'  BK006945.2 NC_001144 1078177   chrXII
13      'Chromosome XIII' BK006946.2 NC_001145  924431   chrXIII
14      'Chromosome XIV'  BK006947.3 NC_001146  784333   chrXIV
15      'Chromosome XV'   BK006948.2 NC_001147 1091291   chrXV
16      'Chromosome XVI'  BK006949.2 NC_001148  948066   chrXVI
17      'Chromosome Mito' AJ011856.1 NC_001224   85779   chrmt
")
names <- c("Chromosome.X", "RefSeq_ID", "Genbank_ID", "chrX")
for (n in names) {
    chr.ref.NC[, n] <- sapply(chr.ref.NC[, n], function(x) as.character(x))
}

#'
#' Return the ref|NC_XXXXXX| of any chromosome number
#'
#' @param chr.num the integer corresponding to the Chromosome
#' @return string of type ref|NC_XXXXXX|
#'
#' @examples
#'
#' GetRefChr(2)
#' # should return 'ref|NC_001134|'
#'
#' @export
GetRefChr <- function(chr.num){
    nc_id <- chr.ref.NC[chr.num, "RefSeq_ID"]
    ref <- paste("ref|", nc_id, "|", sep = "")

    ref
}

#'
#' Return the Length of the chromosome number
#'
#' @param chr.num the integer corresponding to the Chromosome
#'
#' @export
GetLenghtChr <- function(chr.num) {
    chr.ref.NC[chr.ref.NC$chr_num == 2, "Length"]
}

#'
#' Return the chrX of any chromosome reference id (ref|NC_XXXXXX|)
#'
#' @param ref the string corresponding to Refernce id: ref|NC_XXXXXX|
#' @return string of type chrX
#'
#' @export
Ref2Str <- function(ref){
    if (! is.character(ref)) {
        ref <- as.character(ref)
    }
    list.refseq.id <- strsplit(ref, "|", fixed = TRUE)
    refseq.id <- list.refseq.id[[1]][2]
    chrX <- chr.ref.NC[chr.ref.NC$RefSeq_ID == refseq.id, "chrX"]

    chrX
}

#'
#' Return the boolean of Subtelo region
#'
#' @param pos the position (integer)
#' @param chrX the string of the chromosome
#' @param length.subtelo by default 30000
#'
#' @export
Is.InSubtelo <- function(pos, chrX, length.subtelo = 30000) {
    max.length <- chr.ref.NC[chr.ref.NC$chrX == chrX, "Length"]
    # Weird type system of R !!
    pos <- as.integer(pos)
    length.subtelo <- as.integer(length.subtelo)
    if (pos < length.subtelo | pos > (max.length - length.subtelo)) {
        resp <- TRUE
    } else {
        resp <- FALSE
    }

    resp
}

#'
#' Return the Data.Frame updated if the position belong to the Subtelo region
#'
#' @param df the Data.Frame returned by GetWeight
#' @return the data frame df updated
#'
#' @export
GetSubtelo <- function(df) {
    apply(df, 1, function(x) Is.InSubtelo(x[2], x[5]))
}

#'
#' Return the Data.Frame with Random Chromosome and Position
#' (all the weights are set to 1.)
#'
#' @param nsample the number of samples
#' @return the data frame
#'
#' @examples
#'
#' GetRand(10)
#'
#' @export
GetRand <- function(nsample) {
    nchr <- nrow(chr.ref.NC)
    samp <- sample(1:nchr, nsample, replace=TRUE)

    chr <- sapply(samp, function(x) GetRefChr(x))
    orient <- sample(c(-1, 1), nsample, replace=TRUE)

    chrX <- sapply(chr, Ref2Str)
    pos <- sapply(samp, function(x) sample(1:GetLenghtChr(x), 1))

    df <- data.frame(chr = chr
                   , pos = pos
                   , orient = orient
                   , weight = rep(1, nsample)
                   , chrX = chrX
                   , From = rep("Rand", nsample))
                   #, stringAsFactors=TRUE)

    df
}

#'
#' Return
#'
#' @import doBy
#'
#' @export
Is.ValidRand <- function(df.rand) {
    is.ok <- TRUE
    if (max(summaryBy(weight ~ chr + pos + orient, FUN=length, df.rand)[, "weight.length"]) > 1) {
        print("Warning: Same Random position (equal chr, pos, orient).")
        print("         Better to re-launch the random generation.")

        is.ok <- FALSE
    }
    is.ok
}

#'
#' Return
#'
#' @export
GetDataBinded <- function(df) {
    nsample <- nrow(df)
    ra <- GetRand(nsample)
    dat <- rbind(df, ra)

    dat <- data.frame(dat
                    , Subtelo=GetSubtelo(dat))
    dat$From <- as.factor(dat$From)

    dat
}

#'
#' Return
#'
#' @export
GetData <- function(df) {
    df <- GetDataBinded(df)
    df[sample(1:nrow(df), nrow(df)),]
}


#'
#' Return
#'
#' @import pROC
#'
#' @export
GetAUC <- function(df) {
    model <- glm(From ~ Subtelo, weights=weight, df, family="binomial")
    pred <- predict(model, df, type="response")
    roc <- roc(df$From, pred, plot=FALSE)

    auc(roc)
}
