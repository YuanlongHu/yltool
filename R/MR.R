#' run TwoSampleMR
#'
#'
#' @title runTwoSampleMR2
#' @param expID IDs
#' @param outID Outcome IDs
#' @param p1 Significance threshold.
#' @param r2 Clumping r2 cut off.
#' @param kb Clumping distance cutoff.
#' @param maf_threshold MAF threshold to try to infer palindromic SNPs.
#' @param method_list List of methods to use in analysis.
#' @param proxies Look for LD tags?
#' @importFrom TwoSampleMR extract_instruments
#' @importFrom TwoSampleMR extract_outcome_data
#' @importFrom TwoSampleMR harmonise_data
#' @importFrom TwoSampleMR mr
#' @importFrom TwoSampleMR mr_method_list
#' @importFrom pbapply pblapply
#' @return a data frame object
#' @export
#' @author Yuanlong Hu

runTwoSampleMR2 <- function(expID, outID,
                           p1=5e-8, r2=0.01, kb=1000,
                           maf_threshold=0.01,
                           method_list = subset(mr_method_list(), use_by_default)$obj,
                           proxies = TRUE){

  message("*** Find instruments ***")
  data_exp <- pbapply::pblapply(as.list(expID), function(x){
    try({

      exp <- extract_instruments(outcomes=x,
                                 p1=p1, clump=TRUE, r2=r2,
                                 kb=kb, access_token = NULL
      )

      if(is.null(exp)){
        return(data.frame(NULL))
      }else{
        return(exp)
      }
    })
  })


  message("*** Harmonise ***")
  data_exp_out <- pbapply::pblapply(data_exp, function(x){
    try({

      out <- extract_outcome_data(
        snps=x$SNP,
        outcomes=outID,
        proxies = proxies,
        maf_threshold = maf_threshold,
        access_token = NULL
      )

      if(is.null(out)){
        return(data.frame(NULL))
      }else{

        mydata <- harmonise_data(
          exposure_dat=x,
          outcome_dat=out,
          action= 2
        )
        return(mydata)
      }
    })
  })

  message("*** Run MR ***")
  res_exp_out <- pbapply::pblapply(data_exp_out, function(x){

    try({
      res <- mr(x, method_list=method_list)
      return(res)
    })
  })

  res_exp_out <- Reduce(rbind, res_exp_out)
  res_exp_out$pval <- as.numeric(res_exp_out$pval)
  res_exp_out$b <- as.numeric(res_exp_out$b)
  res_exp_out$OR <- exp(res_exp_out$b)

  return(res_exp_out)
}


#' run TwoSampleMR (single)
#'
#'
#' @title runTwoSampleMR
#' @param expID IDs.
#' @param outID Outcome IDs.
#' @param output output type, harmonise or mr.
#' @param p1 Significance threshold.
#' @param r2 Clumping r2 cut off.
#' @param kb Clumping distance cutoff.
#' @param maf_threshold MAF threshold to try to infer palindromic SNPs.
#' @param method_list List of methods to use in analysis.
#' @param proxies Look for LD tags?
#' @importFrom TwoSampleMR extract_instruments
#' @importFrom TwoSampleMR extract_outcome_data
#' @importFrom TwoSampleMR harmonise_data
#' @importFrom TwoSampleMR mr
#' @importFrom TwoSampleMR mr_method_list
#' @importFrom pbapply pblapply
#' @return a MR object
#' @export
#' @author Yuanlong Hu
runTwoSampleMR <- function(expID, outID,output="harmonise",
                            p1=5e-8, r2=0.01, kb=1000,
                            maf_threshold=0.01,
                            method_list = subset(mr_method_list(), use_by_default)$obj,
                            proxies = TRUE){

  message("*** Find instruments ***")

  exp <- extract_instruments(outcomes=expID,
                            p1=p1, clump=T, r2=r2,
                            kb=kb, access_token = NULL)

  message("*** Harmonise ***")
  out <- extract_outcome_data(
        snps=x$SNP,
        outcomes=outID,
        proxies = proxies,
        maf_threshold = maf_threshold,
        access_token = NULL
      )
  mydata <- harmonise_data(
          exposure_dat=x,
          outcome_dat=out,
          action= 2
        )

  message("*** Run MR ***")
  res <- mr(x, method_list=method_list)

  if(output == "harmonise") return(harmonise)
  if(output == "mr") return(res)
}
