\name{NEWS}
\title{trend News}
\encoding{UTF-8}
\section{CHANGES IN TREND VERSION 1.1.4 (2020-09-17)}{
  \subsection{OTHER CHANGES}{
    \itemize{
      \item In function \code{pettitt.test} the p-value
      is now calculated as \code{pval <- min(1, 2.0 * exp(( -6.0 * U^2) / (n^3 + n^2))) to ensure that p does not exceed one.
}
    }
  }
}
\section{CHANGES IN TREND VERSION 1.1.3 (2020-08-31)}{
  \subsection{OTHER CHANGES}{
    \itemize{
      \item CRAN check detected undeclared packages
      that were referred to in Rd files.
      The \sQuote{DESCRIPTION} file was updated 
      with \sQuote{Suggests}.
    }
  }
}
\section{CHANGES IN TREND VERSION 1.1.2 (2020-01-11)}{
  \subsection{OTHER CHANGES}{
    \itemize{
      \item CRAN check detected \sQuote{LTO} related issues.
      The \sQuote{init.c} was re-written.
    }
  }
}
\section{CHANGES IN TREND VERSION 1.1.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The references and description 
      in the .Rd file for the function \code{rrod.test} 
      was extended in order to address the 
      original authors Fligner and Pollicello.
      }
   }
}
\section{CHANGES IN TREND VERSION 1.1.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{lanzante.test} for performing
      a non-parametric single-change point-detection test
      using either the \code{wilcox.test} or the
      robust rank order distribution test (\code{rrod.test})
      was added.
      \item The function \code{rrod.test} for performing
      the non-parametric robust rank order distribution test
      was added.
    }
  }
}
\section{CHANGES IN TREND VERSION 1.0.1}{
  \subsection{MODIFICATION}{
    \itemize{
      \item The internal function \code{.mkScore} was slightly
      recoded to enhance speed. Thanks to Martin Schmelzer.
    }  
  }
}
\section{CHANGES IN TREND VERSION 1.0.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{bartels.test} for performing the rank
      version of von Neumann's ratio test was added.
      \item The function \code{br.test} for performing the
      Buishand range test was added.
      \item The function \code{bu.test} for performing the
      Buishand U test was added.
      \item The function \code{snh.test} for performing Alexandersson's
      standard normal homogeinity test was added.
      \item The generic function \code{plot.cptest} was added
      to plot the output of the above functions.
      \item The function \code{pettitt.test} was recoded in order to
      increase speed.
      \item The functions \code{mk.test}, \code{smk.test},
      \code{cmk.test}, \code{partial.mk.test} and \code{sens.slope}
      were recoded
      in order to increase speed. They consistently produce output
      objects of class \code{"htest"}.
      \item The generic function \code{summary.smktest} was added
      to display seasonal output of \code{smk.test}
      \item The function \code{sea.sens.slope} was recoded and
      produces a scalar output.
      \item The function \code{mult.mk.test} was added for multiside
      or multivariate Mann-Kendall Test.
      \item The \code{man} pages were revisited and now contain
      the underlying mathematical functions.
    }
  }
}
\section{CHANGES IN TREND VERSION 0.2.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{ww.test} for conducting Wald-Wolfowitz
      test for independence and stationarity was added to the package.
    }
  }
}
\section{CHANGES IN TREND VERSION 0.1.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{wm.test} for conducting the Wallis and
      Moore phase-frequency test for randomness was added to the
      package.
      \item the function \code{cs.test} for the Cox-Stuart trend test
      was included.
    }
  }
  \subsection{MODIFICATION}{
    \itemize{
      \item The function \code{sens.slope} was modified to enhance speed.
    }
  }
}
\section{CHANGES IN TREND VERSION 0.0.2}{
  \subsection{BUGFIXES}{
    \itemize{
      \item The NAMESPACE file was updated.
    }
  }
}
