\name{FWD}
\alias{FWD}
\title{Performs MC-ANOVA}
\description{
    This function predicts genetic values drawn from the core using SNPs not in the core and those in the core that are not randomly chosen to be QTL.
}
\usage{
FWD(y, X, df = 20, tol = 1e-7, maxIter = 1000, centerImpute = TRUE,
    verbose = TRUE)
  MC_ANOVA(X, X2 = NULL, core, nQTL, nRep = NULL, maxRep = 300, lambda = 1e-8, sampler = rnorm, ...)
}
\arguments{
    \item{X}{
    }
    \item{X2}{
    }
    \item{core}{
    }
    \item{nQTL}{
    }
    \item{nRep}{
    }
    \item{maxRep}{
    }
    \item{lambda}{
    }
    \item{sampler}{
    }
}
\value{
    A matrix with ancestry groups in the rows, MC-ANOVA correlation estimates in the first column, and Monte Carlo error estimates in the second column.
}