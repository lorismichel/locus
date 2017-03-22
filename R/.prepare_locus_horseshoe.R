prepare_list_hyper_horseshoe <- function(list_hyper, Y, d, p, p_star, q,
                                bool_rmvd_x, bool_rmvd_z,
                                names_x, names_y, names_z, verbose) {

  if (is.null(list_hyper)) {

    if (verbose) cat("list_hyper set automatically. \n")

    list_hyper <- auto_set_hyper_(Y, p, p_star, q)

  } else {

    if (!inherits(list_hyper, c("hyper", "out_hyper")))
      stop(paste("The provided list_hyper must be an object of class ``hyper'' ",
                 "or ``out_hyper''. \n",
                 "*** you must either use the function set_hyper to ",
                 "set your own hyperparameters or use list_hyper from a ``vb'' ",
                 "object or set the argument list_hyper to NULL for automatic choice. ***",
                 sep=""))

    if (inherits(list_hyper, "hyper")) {
      p_hyper_match <- length(bool_rmvd_x)
    } else {
      p_hyper_match <- p
    }


    if (list_hyper$d_hyper != d)
      stop(paste("The dimensions of the provided hyperparameters ",
                 "(list_hyper) are not consistent with that of Y.", sep=""))

    if (list_hyper$p_hyper != p_hyper_match)
      stop(paste("The dimensions of the provided hyperparameters ",
                 "(list_hyper) are not consistent with that of X.", sep=""))

    if (inherits(list_hyper, "hyper")) {
      # remove the entries corresponding to the removed constant predictors in X
      # (if any)
      list_hyper$a <- list_hyper$a[!bool_rmvd_x]
      list_hyper$b <- list_hyper$b[!bool_rmvd_x]
    }

    if (!is.null(names(list_hyper$a)) && names(list_hyper$a) != names_x)
      stop("Provided names for the entries of a do not match the colnames of X.")

    if (!is.null(names(list_hyper$b)) && names(list_hyper$b) != names_x)
      stop("Provided names for the entries of b do not match the colnames of X.")

    if (!is.null(names(list_hyper$eta)) && names(list_hyper$eta) != names_y)
      stop("Provided names for the entries of eta do not match the colnames of Y.")

    if (!is.null(names(list_hyper$kappa)) && names(list_hyper$kappa) != names_y)
      stop("Provided names for the entries of kappa do not match the colnames of Y.")



    if (!is.null(q)) {

      if (inherits(list_hyper, "hyper")) {
        q_hyper_match <- length(bool_rmvd_z)
        # remove the entries corresponding to the removed constant predictors in X
        # (if any)
        list_hyper$phi <- list_hyper$phi[!bool_rmvd_z]
        list_hyper$xi <- list_hyper$xi[!bool_rmvd_z]
      } else {
        q_hyper_match <- q
      }

      if (list_hyper$q_hyper != q_hyper_match)
        stop(paste("The dimensions of the provided hyperparameters ",
                   "(list_hyper) are not consistent with that of Z.", sep=""))

      if (!is.null(names(list_hyper$phi)) && names(list_hyper$phi) != names_z)
        stop("Provided names for the entries of phi do not match the colnames of Z.")

      if (!is.null(names(list_hyper$xi)) && names(list_hyper$xi) != names_z)
        stop("Provided names for the entries of xi do not match the colnames of Z.")

    }

  }

  class(list_hyper) <- "out_hyper"

  list_hyper
}


prepare_list_init_ <- function(list_init, Y, d, p, p_star, q, bool_rmvd_x,
                               bool_rmvd_z, names_x, names_y, names_z,
                               user_seed, verbose) {

  if (is.null(list_init)) {

    if (!is.null(user_seed) & verbose) cat(paste("Seed set to user_seed ",
                                                 user_seed,". \n", sep=""))

    if (verbose) cat(paste("list_init set automatically. \n", sep=""))

    list_init <- auto_set_init_(Y, p, p_star, user_seed, q)

  } else {

    if (!is.null(user_seed))
      warning("user_seed not used since a non-NULL list_init was provided. \n")

    if (!inherits(list_init, c("init", "out_init")))
      stop(paste("The provided list_init must be an object of class ``init'' or ",
                 " `` out_init''. \n",
                 "*** you must either use the function set_init to ",
                 "set your own initialization or use list_init from a ``vb'' ",
                 "object or  set the argument list_init to NULL for automatic ",
                 "initialization. ***",
                 sep=""))

    if (inherits(list_init, "init")) {
      p_init_match <- length(bool_rmvd_x)
    } else {
      p_init_match <- p
    }


    if (list_init$d_init != d)
      stop(paste("The dimensions of the provided initial parameters ",
                 "(list_init) are not consistent with that of Y.", sep=""))

    if (list_init$p_init != p_init_match)
      stop(paste("The dimensions of the provided initial parameters ",
                 "(list_init) are not consistent with that of X.", sep=""))

    if (inherits(list_init, "init")) {
      # remove the entries corresponding to the removed constant predictors in X
      # (if any)
      list_init$b_vb <- list_init$b_vb[!bool_rmvd_x,, drop = FALSE]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[!bool_rmvd_x,, drop = FALSE]
    }

    if (!is.null(q)) {

      if (inherits(list_init, "init")) {
        q_init_match <- length(bool_rmvd_z)
        # remove the entries corresponding to the removed constant predictors in X
        # (if any)
        list_init$mu_alpha_vb <- list_init$mu_alpha_vb[!bool_rmvd_z,, drop = FALSE]
        list_init$sig2_alpha_vb <- list_init$sig2_alpha_vb[!bool_rmvd_z,, drop = FALSE]
      } else {
        q_init_match <- q
      }

      if (list_init$q_init != q_init_match)
        stop(paste("The dimensions of the provided initial parameters ",
                   "(list_init) are not consistent with that of Z.", sep=""))
    }

  }

  class(list_init) <- "out_init"

  list_init
}
