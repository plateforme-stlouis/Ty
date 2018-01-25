
# from doc: http://samtools.github.io/hts-specs/SAMv1.pdf

fake.bit <- 2 ** (0:11)
from.doc <- c(
    "template having multiple segments in sequencing"
  , "each segment properly aligned according to the aligner"
  , "segment unmapped"
  , "next segment in the template unmapped"
  , "SEQ being reverse complemented"
  , "SEQ of the next segment in the template being reverse complemented"
  , "the first segment in the template"
  , "the last segment in the template"
  , "secondary alignment"
  , "not passing filters, such as platform/vendor quality controls"
  , "PCR or optical duplicate"
  , "supplementary alignment"
)
description <- c(
    "Read paired"
  , "Read mapped in proper pair"
  , "Read unmapped"
  , "Mate unmapped"
  , " Read reverse strand"
  , "Mate reverse strand"
  , "First in pair"
  , "Second in pair"
  , "Not primary alignment"
  , " Read fails"
  , "Read is PCR or optical duplicate"
  , "Supplementary alignment"
)
df.flag <- data.frame(fake.bit, description, stringsAsFactors = FALSE)

GetFakeBit <- function(flag) {
    fb <- df.flag$fake.bit
    that <- NULL

    while (length(fb) > 0) {
        fb <- fb[fb <= flag]
        m <- fb[length(fb)]

        that <- c(that, m)
        flag <- flag - m
    }

    that
}

InfoFlag <- function(flag) {
    subset(df.flag, df.flag$fake.bit %in% GetFakeBit(flag))
}


Swap <- function(lst, bit1, bit2){
    is.bit1.in <- bit1 %in% lst
    is.bit2.in <- bit2 %in% lst
    if (is.bit1.in && is.bit2.in) {
        i <- lst == bit1
        j <- lst == bit2
        lst[i] <- bit2
        lst[j] <- bit1
    } else {
        if (is.bit1.in) {
            lst[lst == bit1] <- bit2
        }
        if (is.bit2.in) {
            lst[lst == bit2] <- bit1
        }
    }
    lst
}

RevFakeBit <- function(flag) {
    that <- GetFakeBit(flag)

    that <- Swap(that, 4, 8)
    that <- Swap(that, 16, 32)
    that <- Swap(that, 64, 128)

    that
}

RevFlag <- function(flag) {
    Reduce(sum, RevFakeBit(flag))
}

InfoRevFlag <- function(flag) {
    InfoFlag(RevFlag(flag))
}




FlagToBits <- function(int) {
    b <- intToBits(int)
    n <- length(b)
    res <- rep(0, n)
    for (i in 1:n){
        if (b[i] == 1) {
            res[i] <- 1
        }
    }
    res
}

FlagInfo <- function(flag) {
    subset(df.flag, FlagToBits(flag) == 1)
}


BitToFlag <- function(b) {
    f <- 0
    for (i in 1:length(b)) {
        if (b[i] == 1) {
            f <- f + 2 ** (i - 1)
        }
    }
    f
}

SwapByPos <- function(lst, pos1, pos2) {
    pos1 <- pos1 + 1
    pos2 <- pos2 + 1
    is.pos1.in <- lst[pos1] == 1
    is.pos2.in <- lst[pos2] == 1
    if (is.pos1.in && is.pos2.in) {
        return(lst)
    } else {
        if (is.pos1.in) {
            lst[pos1] <- 0
            lst[pos2] <- 1
        }
        if (is.pos2.in) {
            lst[pos1] <- 1
            lst[pos2] <- 0
        }
    }
    lst
}

RevBit <- function(flag){
    that <- FlagToBits(flag)

    that <- SwapByPos(that, 2, 3)
    that <- SwapByPos(that, 4, 5)
    that <- SwapByPos(that, 6, 7)

    that
}

FlagRev <- function(flag) {
    BitToFlag(RevBit(flag))

}

FlagRevInfo <- function(flag) {
    FlagInfo(FlagRev(flag))
}
