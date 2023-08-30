library(edgeR)

# Useful for situations where you need to batch correct a dataset where
# there is partial confounding between treatments and batches, but there
# are some conditions, eg controls, that are in all batches.
# Feed the shared conditions into the counts, batch, and covar_mod argument
# and the full dataset into batch_all, counts_all, covar_mod_all
# Batch effect is learned from the shared conditions and applied to the full dataset

combat_with_reference <- function (counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE, 
    shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL,
    batch_all = NULL, counts_all = NULL, covar_mod_all = NULL
) {
    batch <- as.factor(batch)
    if (any(table(batch) <= 1)) {
        stop("ComBat-seq doesn't support 1 sample per batch yet")
    }
    keep_lst <- lapply(levels(batch), function(b) {
        which(apply(counts[, batch == b], 1, function(x) {
            !all(x == 0)
        }))
    })
    keep <- Reduce(intersect, keep_lst)
    rm <- setdiff(1:nrow(counts), keep)
    countsOri <- counts
    counts <- counts[keep, ]
    dge_obj <- DGEList(counts = counts)
    n_batch <- nlevels(batch)
    batches_ind <- lapply(1:n_batch, function(i) {
        which(batch == levels(batch)[i])
    })
    n_batches <- sapply(batches_ind, length)
    n_sample <- sum(n_batches)
    batches_ind_all <- lapply(1:n_batch, function(i) {
        which(batch_all == levels(batch)[i])
    })
    n_batches_all <- sapply(batches_ind_all, length)
    n_sample_all <- sum(n_batches_all)
    cat("Found", n_batch, "batches\n")
    batchmod <- model.matrix(~-1 + batch)
    group <- as.factor(group)
    if (full_mod & nlevels(group) > 1) {
        cat("Using full model in ComBat-seq.\n")
        mod <- model.matrix(~group)
    }
    else {
        cat("Using null model in ComBat-seq.\n")
        mod <- model.matrix(~1, data = as.data.frame(t(counts)))
    }
    if (!is.null(covar_mod)) {
        if (is.data.frame(covar_mod)) {
            covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), 
                function(i) {
                  model.matrix(~covar_mod[, i])
                }))
        }
        covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x) {
            all(x == 1)
        })]
    }
    mod <- cbind(mod, covar_mod)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
    if (qr(design)$rank < ncol(design)) {
        if (ncol(design) == (n_batch + 1)) {
            stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")
        }
        if (ncol(design) > (n_batch + 1)) {
            if ((qr(design[, -c(1:n_batch)])$rank < ncol(design[, 
                -c(1:n_batch)]))) {
                stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
            }
            else {
                stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")
            }
        }
    }
    NAs = any(is.na(counts))
    if (NAs) {
        cat(c("Found", sum(is.na(counts)), "Missing Data Values\n"), 
            sep = " ")
    }
    cat("Estimating dispersions\n")
    disp_common <- sapply(1:n_batch, function(i) {
        if ((n_batches[i] <= ncol(design) - ncol(batchmod) + 
            1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)) {
            return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                design = NULL, subset = nrow(counts)))
        }
        else {
            return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], 
                design = mod[batches_ind[[i]], ], subset = nrow(counts)))
        }
    })
    genewise_disp_lst <- lapply(1:n_batch, function(j) {
        if ((n_batches[j] <= ncol(design) - ncol(batchmod) + 
            1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)) {
            return(rep(disp_common[j], nrow(counts)))
        }
        else {
            return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], 
                design = mod[batches_ind[[j]], ], dispersion = disp_common[j], 
                prior.df = 0))
        }
    })
    names(genewise_disp_lst) <- paste0("batch", levels(batch))
    phi_matrix <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
    for (k in 1:n_batch) {
        phi_matrix[, batches_ind[[k]]] <- sva:::vec2mat(genewise_disp_lst[[k]], 
            n_batches[k])
    }
    cat("Fitting the GLM model\n")
    glm_f <- glmFit(dge_obj, design = design, dispersion = phi_matrix, 
        prior.count = 1e-04)
    alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample)
    new_offset <- t(sva:::vec2mat(getOffset(dge_obj), nrow(counts))) + 
        sva:::vec2mat(alpha_g, ncol(counts))
    glm_f2 <- glmFit.default(dge_obj$counts, design = design, 
        dispersion = phi_matrix, offset = new_offset, prior.count = 1e-04)
    gamma_hat <- glm_f2$coefficients[, 1:n_batch]
    mu_hat <- glm_f2$fitted.values
    phi_hat <- do.call(cbind, genewise_disp_lst)
    counts_all_ori <- counts_all
    counts_all <- counts_all[keep, ]
    dge_obj_all <- DGEList(counts = counts_all)
    new_offset_all <- t(sva:::vec2mat(getOffset(dge_obj_all), nrow(counts_all))) + 
        sva:::vec2mat(alpha_g, ncol(counts_all))
    phi_matrix_all <- matrix(NA, nrow = nrow(counts_all), ncol = ncol(counts_all))
    for (k in 1:n_batch) {
        phi_matrix_all[, batches_ind_all[[k]]] <- sva:::vec2mat(genewise_disp_lst[[k]], 
            n_batches_all[k])
    }
    glm_all <- glmFit.default(
        counts_all, design = covar_mod_all, dispersion = phi_matrix_all, 
        offset = new_offset_all, prior.count = 1e-04
    )
    mu_hat_all <- glm_all$fitted.values
    if (shrink) {
        stop("Not implemented")
    }
    else {
        cat("Shrinkage off - using GLM estimates for parameters\n")
        gamma_star_mat <- gamma_hat
        phi_star_mat <- phi_hat
    }
    mu_star <- matrix(NA, nrow = nrow(counts_all), ncol = ncol(counts_all))
    for (jj in 1:n_batch) {
        mu_star[, batches_ind_all[[jj]]] <- exp(log( mu_hat_all[, batches_ind_all[[jj]]]) - 
            sva:::vec2mat(gamma_star_mat[, jj], n_batches_all[jj]))
    }
    phi_star <- rowMeans(phi_star_mat)
    cat("Adjusting the data\n")
    adjust_counts <- matrix(NA, nrow = nrow(counts_all), ncol = ncol(counts_all))
    for (kk in 1:n_batch) {
        counts_sub <- counts_all[, batches_ind_all[[kk]]]
        old_mu <- mu_hat_all[, batches_ind_all[[kk]]]
        old_phi <- phi_hat[, kk]
        new_mu <- mu_star[, batches_ind_all[[kk]]]
        new_phi <- phi_star
        adjust_counts[, batches_ind_all[[kk]]] <- sva:::match_quantiles(
            counts_sub = counts_sub, 
            old_mu = old_mu, old_phi = old_phi, new_mu = new_mu, 
            new_phi = new_phi)
    }
    adjust_counts_whole <- matrix(NA, nrow = nrow(counts_all_ori), 
        ncol = ncol(counts_all_ori))
    dimnames(adjust_counts_whole) <- dimnames(counts_all_ori)
    adjust_counts_whole[keep, ] <- adjust_counts
    adjust_counts_whole[rm, ] <- counts_all_ori[rm, ]
    return(adjust_counts_whole)
}
