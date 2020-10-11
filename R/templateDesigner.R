#' Check sequence
#'
#' This function inputs a character string and checks that all values
#' are ATCG
#'
#' @param char.input A DNA sequence as a string containing ATCG
#' @return A boolean value stating whether the sequence is valid
#' @export
validate.sequence <- function(char.input) {
  if (stringi::stri_detect_charclass(toupper(char.input), "[^ATCG]", negate = TRUE)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Remove whitespace from DNA sequence
#'
#' This function removes whitespace from the sequence string
#'
#' @param seq A DNA sequence as a string containing ATCG
#' @return TRUE/FALSE stating whether the sequence is valid
#' @export
clean.input <- function(seq) {
  clean.seq <- stringi::stri_replace_all_charclass(seq, "\\p{WHITE_SPACE}", "")
}

#' Find PAM sites and create sgRNA
#'
#' This scans the DNA sequence for PAM sites and builds sgRNA guides
#' for the sense and antisense strand
#'
#' @param seq A DNA sequence as a string containing ATCG
#' @return a dataframe with the direction of the guide, index and sequence
#' @export
find.pams <- function(seq) {
  seq <- clean.input(seq)
  if (nchar(seq) < 150 | !validate.sequence(seq)) {
    stop("Please input a valid sequence at least 150bp")
  } else{
    pam.res <- tryCatch({
      res <- NULL
      full.seq <- seq
      sub.seq <- substr(full.seq, 60, nchar(full.seq) - 60)
      forward.pam <-
        unlist(gregexpr("(?=[A-z]GG)", toupper(sub.seq), perl = TRUE))
      if (!forward.pam[1] == -1L) {
        for (i in 1:length(forward.pam)) {
          stop <- forward.pam[i] + 61
          sgRNA.f <- substr(full.seq, stop - 22, stop)
          res <-
            rbind(res,
                  data.frame(
                    dir = "F",
                    index = stop,
                    guide = sgRNA.f
                  ))
        }
      }
      reverse.pam <-
        unlist(gregexpr("(?=CC[A-z])", toupper(sub.seq), perl = TRUE))
      if (!reverse.pam[1] == -1L)
        for (j in 1:length(reverse.pam)) {
          start <- reverse.pam[j] + 59
          sgRNA.r <- substr(full.seq, start, start + 22)
          res <-
            rbind(res,
                  data.frame(
                    dir = "R",
                    index = start,
                    guide = sgRNA.r
                  ))
        }
      return(res)
    },
    error = function(cond) {
      return(FALSE)
    })
    return(pam.res)
  }
}

#' Create a donor oligo based on input sequence and sgRNA guides
#'
#' This function finds 60bp homology arms up and downstream of the sgRNA cut site and
#' merges the left homology arm, insert and right homology arm to build a donor template
#'
#' @param pam.res The dataframe results after running find.pams()
#' @param full.seq A sequence 150bp or longer
#' @param ins A DNA sequence that acts as the insert
#' @return a dataframe with the guide sequence and the corresponding donor oligo
#' @export
get.homology.template <- function(full.seq, ins) {
  full.seq <- clean.input(full.seq)
  ins <- clean.input(ins)
  pam.res <- find.pams(full.seq)
  homology.arms <- tryCatch({
    hdr.template <- NULL
    for (i in 1:nrow(pam.res)) {
      guide <- pam.res$guide[i]
      index <- pam.res$index[i]
      if (pam.res$dir[i] == "F") {
        left.arm <- substring(full.seq,index-66, index - 7)
        right.arm <- substring(full.seq, index - 6, index + 53)
      } else if (pam.res$dir[i] == "R") {
        left.arm <- substring(full.seq,  index-53 , index + 6)
        right.arm <- substring(full.seq, index + 7, index + 66)
      }
      hdr.template <- rbind(hdr.template,data.frame(guide,paste0(left.arm, ins , right.arm)))

    }
    colnames(hdr.template) <- c("sgRNA Guide", "Oligo donor")
    return(hdr.template)
  },
  error = function(cond) {
    return(FALSE)
  })
}
