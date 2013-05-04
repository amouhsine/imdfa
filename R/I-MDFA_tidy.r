#library(formatR)
#tidy.source(source="imdfa/inst/I-MDFA_new14.r", keep.comment=FALSE, keep.blank.line=FALSE)
spec_mat_comp <- function(weight_func, L, Lag) {
    K <- length(weight_func[, 1]) - 1
    weight_h <- weight_func
    weight_h[1, ] <- weight_h[1, ] * 0.5
    weight_target <- weight_h[, 1]
    weight_h <- weight_h * exp(-(0+1i) * Arg(weight_target))
    weight_target <- weight_target * exp(-(0+1i) * Arg(weight_target))
    weight_h_exp <- as.matrix(weight_h[, 2:(dim(weight_h)[2])])
    spec_mat <- as.vector(t(as.matrix(weight_h_exp[1, ]) %*% t(as.matrix(rep(1, L)))))
    for (j in 1:(K)) {
        omegak <- j * pi/K
        exp_vec <- exp((0+1i) * omegak * ((0:(L - 1)) - Lag))
        spec_mat <- cbind(spec_mat, as.vector(t(as.matrix(weight_h_exp[j + 1, ]) %*% t(as.matrix(exp_vec)))))
    }
    dim(spec_mat)
    return(list(spec_mat = spec_mat))
}
mat_func <- function(i1, i2, L, weight_h_exp, lambda_decay, lambda_cross, lambda_smooth, 
        Lag, weight_constraint, shift_constraint, grand_mean) {
    if (Lag > (L - 1)/2) {
        print("Lag larger than L/2!!!!! Will be trimmed automtically to L/2 (symmetric filter)")
        Lag <- as.integer(L/2)
    }
    Q_smooth <- matrix(data = rep(0, ((L) * length(weight_h_exp[1, ]))^2), 
            nrow = (L) * length(weight_h_exp[1, ]), 
            ncol = (L) * length(weight_h_exp[1, ]))
    Q_decay <- matrix(data = rep(0, ((L) * length(weight_h_exp[1, ]))^2), 
            nrow = (L) * length(weight_h_exp[1, ]), 
            ncol = (L) * length(weight_h_exp[1, ]))
    if ((length(weight_h_exp[1, ]) > 1)) {
        Q_cross <- matrix(data = rep(0, ((L) * length(weight_h_exp[1, ]))^2), 
                nrow = (L) * length(weight_h_exp[1, ]), 
                ncol = (L) * length(weight_h_exp[1, ]))
        Q_centraldev_original <- matrix(data = rep(0, ((L) * length(weight_h_exp[1,]))^2), nrow = (L) * length(weight_h_exp[1, ]), ncol = (L) * length(weight_h_exp[1,]))
    }
    else {
        Q_cross <- NULL
    }
    for (i in 1:L) {
        Q_decay[i, i] <- (1 + lambda_decay[1])^(2 * abs(i - 1 - max(0, Lag)))
        if (i == 1) {
            Q_smooth[i, i:(i + 2)] <- c(1, -2, 1)
        }
        else {
            if (i == 2) {
                Q_smooth[i, (i - 1):(i + 2)] <- c(-2, 5, -4, 1)
            }
            else {
                if (i == L) {
                    Q_smooth[i, (i - 2):i] <- c(1, -2, 1)
                }
                else {
                    if (i == L - 1) {
                        Q_smooth[i, (i - 2):(i + 1)] <- c(1, -4, 5, -2)
                    }
                    else {
                        Q_smooth[i, (i - 2):(i + 2)] <- c(1, -4, 6, -4, 1)
                    }
                }
            }
        }
    }
    if (length(weight_h_exp[1, ]) > 1) {
        for (j in 1:max(1, (length(weight_h_exp[1, ]) - 1))) {
            Q_smooth[j * L + 1:L, j * L + 1:L] <- Q_smooth[1:L, 1:L]
            Q_decay[j * L + 1:L, j * L + 1:L] <- Q_decay[1:L, 1:L]
        }
        Q_centraldev_original <- diag(rep(1, L * length(weight_h_exp[1, ])))
        diag(Q_centraldev_original[1:L, L + 1:L]) <- rep(-1, L)
        for (i in 2:length(weight_h_exp[1, ])) {
            diag(Q_centraldev_original[(i - 1) * L + 1:L, 1:L]) <- rep(1, L)
            diag(Q_centraldev_original[(i - 1) * L + 1:L, (i - 1) * L + 1:L]) <- rep(1, L)
            diag(Q_centraldev_original[1:L, (i - 1) * L + 1:L]) <- rep(-1, L)
        }
        Q_centraldev_original <- solve(Q_centraldev_original)
        if (grand_mean) {
            diag(Q_cross[L + 1:((length(weight_h_exp[1, ]) - 1) * L), L + 1:((length(weight_h_exp[1,]) - 1) * L)]) <- rep(1, ((length(weight_h_exp[1, ]) - 1) * L))
        }
        else {
            diag(Q_cross) <- 1
            for (i in 1:length(weight_h_exp[1, ])) {
                for (j in 1:L) {
                    Q_cross[(i - 1) * L + j, j + (0:(length(weight_h_exp[1, ]) - 1)) * 
                                    L] <- Q_cross[(i - 1) * L + j, j + (0:(length(weight_h_exp[1,]) - 1)) * L] - 1/length(weight_h_exp[1, ])
                }
            }
        }
    }
    else {
        Q_centraldev_original <- NULL
    }
    Q_decay <- Q_decay * lambda_decay[2]
    Q_cross <- Q_cross * lambda_cross
    Q_smooth <- Q_smooth * lambda_smooth
    if (lambda_decay[2] > 0) {
        Q_decay <- lambda_decay[2] * (Q_decay/(sum(diag(Q_decay))))
    }
    if (lambda_cross > 0) {
        Q_cross <- lambda_cross * (Q_cross/(sum(diag(Q_cross))))
    }
    if (lambda_smooth > 0) {
        Q_smooth <- lambda_smooth * (Q_smooth/(sum(diag(Q_smooth))))
    }
    if (i1) {
        if (i2) {
            if (Lag < 1) {
                w_eight <- c(-(Lag - 1) * weight_constraint[1] - shift_constraint[1], 
                        Lag * weight_constraint[1] + shift_constraint[1], rep(0, L - 2))
            }
            else {
                w_eight <- c(rep(0, Lag), weight_constraint[1] - shift_constraint[1], 
                        shift_constraint[1], rep(0, L - Lag - 2))
            }
            if (length(weight_h_exp[1, ]) > 1) {
                for (j in 2:length(weight_h_exp[1, ])) {
                    if (Lag < 1) {
                        w_eight <- c(w_eight, -(Lag - 1) * weight_constraint[j] - shift_constraint[j], 
                                Lag * weight_constraint[j] + shift_constraint[j], rep(0, L - 2))
                    }
                    else {
                        w_eight <- c(w_eight, c(rep(0, Lag), weight_constraint[j] - shift_constraint[j], 
                                        shift_constraint[j], rep(0, L - Lag - 2)))
                    }
                }
            }
        }
        else {
            if (Lag < 1) {
                w_eight <- c(weight_constraint[1], rep(0, L - 1))
            }
            else {
                w_eight <- c(rep(0, Lag), weight_constraint[1], rep(0, L - Lag - 
                                        1))
            }
            if (length(weight_h_exp[1, ]) > 1) {
                for (j in 2:length(weight_h_exp[1, ])) {
                    if (Lag < 1) {
                        w_eight <- c(w_eight, weight_constraint[j], rep(0, L - 1))
                    }
                    else {
                        w_eight <- c(w_eight, rep(0, Lag), weight_constraint[j], rep(0, L - Lag - 1))
                    }
                }
            }
        }
    }
    else {
        if (i2) {
            if (Lag < 1) {
                w_eight <- c(0, shift_constraint[1]/(1 - Lag), rep(0, L - 2))
            }
            else {
                w_eight <- c(rep(0, Lag + 1), shift_constraint[1], rep(0, L - Lag - 2))
            }
            if (length(weight_h_exp[1, ]) > 1) {
                for (j in 2:length(weight_h_exp[1, ])) {
                    if (Lag < 1) {
                        w_eight <- c(w_eight, c(0, shift_constraint[j]/(1 - Lag), rep(0, L - 2)))
                    }
                    else {
                        w_eight <- c(w_eight, c(rep(0, Lag + 1), shift_constraint[j], rep(0, L - Lag - 2)))
                    }
                }
            }
        }
        else {
            w_eight <- rep(0, L * length(weight_h_exp[1, ]))
        }
    }
    if (i2) {
        if (i1) {
            des_mat <- matrix(data = rep(0, (L - 2) * L * (length(weight_h_exp[1,]))^2), 
                    nrow = (L - 2) * length(weight_h_exp[1, ]), 
                    ncol = (L) * length(weight_h_exp[1, ]))
            for (i in 1:(L - 2)) {
                if (Lag < 1) {
                    des_mat[i, i + 2 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- 1
                    des_mat[i, 1 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- i
                    des_mat[i, 2 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- -(i + 1)
                }
                else {
                    des_mat[i, ifelse(i < Lag + 1, i, i + 2) + (0:(length(weight_h_exp[1,]) - 1)) * L] <- 1
                    des_mat[i, Lag + 1 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- ifelse(i < Lag + 1, -(Lag + 2 - i), i - Lag)
                    des_mat[i, Lag + 2 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- ifelse(i < Lag + 1, (Lag + 1 - i), -(i - Lag + 1))
                }
            }
            if (length(weight_h_exp[1, ]) > 1) {
                for (j in 1:max(1, (length(weight_h_exp[1, ]) - 1))) {
                    for (i in 1:(L - 2)) {
                        if (Lag < 1) {
                            des_mat[i + j * (L - 2), i + 2] <- -1
                            des_mat[i + j * (L - 2), 1] <- -i
                            des_mat[i + j * (L - 2), 2] <- (i + 1)
                        }
                        else {
                            des_mat[i + j * (L - 2), ifelse(i < Lag + 1, i, i + 2)] <- -1
                            des_mat[i + j * (L - 2), Lag + 1] <- -ifelse(i < Lag + 1, -(Lag + 2 - i), i - Lag)
                            des_mat[i + j * (L - 2), Lag + 2] <- -ifelse(i < Lag + 1, (Lag + 1 - i), -(i - Lag + 1))
                        }
                        if (Lag < 1) {
                            des_mat[i + j * (L - 2), i + 2 + j * L] <- 1
                            des_mat[i + j * (L - 2), 1 + j * L] <- i
                            des_mat[i + j * (L - 2), 2 + j * L] <- -(i + 1)
                        }
                        else {
                            des_mat[i + j * (L - 2), ifelse(i < Lag + 1, i + j * L, i + 2 + j * L)] <- 1
                            des_mat[i + j * (L - 2), Lag + 1 + j * L] <- ifelse(i < Lag + 1, -(Lag + 2 - i), i - Lag)
                            des_mat[i + j * (L - 2), Lag + 2 + j * L] <- ifelse(i < Lag + 1, (Lag + 1 - i), -(i - Lag + 1))
                        }
                    }
                }
            }
        }
        else {
            des_mat <- matrix(data = rep(0, (L - 1) * L * (length(weight_h_exp[1,]))^2), 
                    nrow = (L - 1) * length(weight_h_exp[1, ]), 
                    ncol = (L) * length(weight_h_exp[1, ]))
            for (i in 1:(L - 1)) {
                if (Lag < 1) {
                    des_mat[i, ifelse(i < 2, i, i + 1) + (0:(length(weight_h_exp[1,]) - 1)) * L] <- 1
                    des_mat[i, 1 + 1 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- ifelse(i == 1, Lag, Lag - i)/(1 - Lag)
                }
                else {
                    des_mat[i, ifelse(i < Lag + 2, i, i + 1) + (0:(length(weight_h_exp[1,]) - 1)) * L] <- 1
                    des_mat[i, Lag + 2 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- ifelse(i <  Lag + 2, Lag + 1 - i, Lag - i)
                }
            }
            if (length(weight_h_exp[1, ]) > 1) {
                for (j in 1:max(1, (length(weight_h_exp[1, ]) - 1))) {
                    for (i in 1:(L - 1)) {
                        if (Lag < 1) {
                            des_mat[i + j * (L - 1), ifelse(i < 2, i, i + 1)] <- -1
                            des_mat[i + j * (L - 1), 1 + 1] <- -ifelse(i == 1, Lag, Lag - 
                                            i)/(1 - Lag)
                        }
                        else {
                            des_mat[i + j * (L - 1), ifelse(i < Lag + 2, i, i + 1)] <- -1
                            des_mat[i + j * (L - 1), Lag + 2] <- -ifelse(i < Lag + 2, Lag + 1 - i, Lag - i)
                        }
                        if (Lag < 1) {
                            des_mat[i + j * (L - 1), ifelse(i < 2, i, i + 1) + j * L] <- 1
                            des_mat[i + j * (L - 1), 1 + 1 + j * L] <- ifelse(i == 1, Lag, Lag - i)/(1 - Lag)
                        }
                        else {
                            des_mat[i + j * (L - 1), ifelse(i < Lag + 2, i, i + 1) + j * L] <- 1
                            des_mat[i + j * (L - 1), Lag + 2 + j * L] <- ifelse(i < Lag + 2, Lag + 1 - i, Lag - i)
                        }
                    }
                }
            }
        }
    }
    else {
        if (i1) {
            des_mat <- matrix(data = rep(0, (L - 1) * L * (length(weight_h_exp[1, ]))^2), 
                    nrow = (L - 1) * length(weight_h_exp[1, ]), 
                    ncol = (L) * length(weight_h_exp[1, ]))
            for (i in 1:(L - 1)) {
                if (Lag < 1) {
                    des_mat[i, i + 1 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- 1
                    des_mat[i, 1 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- -1
                }
                else {
                    des_mat[i, ifelse(i < Lag + 1, i, i + 1) + (0:(length(weight_h_exp[1,]) - 1)) * L] <- 1
                    des_mat[i, Lag + 1 + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- -1
                }
            }
            if (length(weight_h_exp[1, ]) > 1) {
                for (j in 1:max(1, (length(weight_h_exp[1, ]) - 1))) {
                    for (i in 1:(L - 1)) {
                        if (Lag < 1) {
                            des_mat[i + j * (L - 1), i + 1] <- -1
                            des_mat[i + j * (L - 1), 1] <- 1
                            des_mat[i + j * (L - 1), i + 1 + j * L] <- 1
                            des_mat[i + j * (L - 1), 1 + j * L] <- -1
                        }
                        else {
                            des_mat[i + j * (L - 1), ifelse(i < Lag + 1, i, i + 1)] <- -1
                            des_mat[i + j * (L - 1), Lag + 1] <- 1
                            des_mat[i + j * (L - 1), ifelse(i < Lag + 1, i, i + 1) + j * L] <- 1
                            des_mat[i + j * (L - 1), Lag + 1 + j * L] <- -1
                        }
                    }
                }
            }
        }
        else {
            des_mat <- matrix(data = rep(0, (L) * L * (length(weight_h_exp[1, ]))^2), 
                    nrow = (L) * length(weight_h_exp[1, ]), 
                    ncol = (L) * length(weight_h_exp[1, ]))
            for (i in 1:(L)) {
                des_mat[i, i + (0:(length(weight_h_exp[1, ]) - 1)) * L] <- 1
            }
            if (length(weight_h_exp[1, ]) > 1) {
                for (j in 1:max(1, (length(weight_h_exp[1, ]) - 1))) {
                    for (i in 1:(L)) {
                        des_mat[i + (j) * (L), i] <- -1
                        des_mat[i + (j) * (L), i + j * L] <- 1
                    }
                }
            }
        }
    }
    if (!grand_mean & length(weight_h_exp[1, ]) > 1) {
        des_mat <- t(Q_centraldev_original %*% t(des_mat))
    }
    if ((length(weight_h_exp[1, ]) > 1)) {
        if (grand_mean) {
            reg_t <- (Q_smooth + Q_decay + t(Q_centraldev_original) %*% Q_cross %*% Q_centraldev_original)
        }
        else {
            reg_t <- (Q_smooth + Q_decay + Q_cross)
        }
    }
    else {
        reg_t <- (Q_smooth + Q_decay)
    }
    reg_mat <- (des_mat) %*% reg_t %*% t(des_mat)
    if (lambda_smooth + lambda_decay[2] + lambda_cross > 0) {
        disentangle_des_mat_effect <- sum(diag(reg_t))/sum(diag(reg_mat))
        reg_mat <- reg_mat * disentangle_des_mat_effect
        reg_xtxy <- des_mat %*% reg_t %*% w_eight * (disentangle_des_mat_effect)
    }
    else {
        reg_xtxy <- des_mat %*% reg_t %*% w_eight
    }
    return(list(des_mat = des_mat, reg_mat = reg_mat, reg_xtxy = reg_xtxy, w_eight = w_eight))
}
mdfa_analytic_new <- function(K, L, lambda, weight_func, Lag, Gamma, expweight, cutoff, 
        i1, i2, weight_constraint, lambda_cross, lambda_decay, lambda_smooth, lin_expweight, 
        shift_constraint, grand_mean) {
    spec_mat <- spec_mat_comp(weight_func, L, Lag)$spec_mat
    omega_Gamma <- as.integer(cutoff * K/pi)
    if ((K - omega_Gamma + 1) > 0) {
        if (lin_expweight) {
            expweight_vec <- c(rep(1, omega_Gamma), 1 + rep(expweight, K - omega_Gamma + 1))
        }
        else {
            expweight_vec <- c(rep(1, omega_Gamma), (1:(K - omega_Gamma + 1))^(expweight/2))
        }
        weight_h <- weight_func * expweight_vec
    }
    else {
        expweight_vec <- rep(1, K + 1)
        weight_h <- weight_func * expweight_vec
    }
    weight_h[1, ] <- weight_h[1, ] * 0.5
    weight_target <- weight_h[, 1]
    weight_h <- weight_h * exp(-(0+1i) * Arg(weight_target))
    weight_target <- Re(weight_target * exp(-(0+1i) * Arg(weight_target)))
    weight_h_exp <- as.matrix(weight_h[, 2:(dim(weight_h)[2])])
    spec_mat <- t(t(spec_mat) * expweight_vec)
    mat_obj <- mat_func(i1, i2, L, weight_h_exp, lambda_decay, lambda_cross, lambda_smooth, 
            Lag, weight_constraint, shift_constraint, grand_mean)
    des_mat <- mat_obj$des_mat
    reg_mat <- mat_obj$reg_mat
    reg_xtxy <- mat_obj$reg_xtxy
    w_eight <- mat_obj$w_eight
    mat_x <- des_mat %*% spec_mat
    X_new <- t(Re(mat_x)) + sqrt(1 + Gamma * lambda) * (0+1i) * t(Im(mat_x))
    xtx <- t(Re(X_new)) %*% Re(X_new) + t(Im(X_new)) %*% Im(X_new)
    xtxy <- t(Re(t(w_eight) %*% spec_mat) %*% Re(t(spec_mat) %*% t(des_mat)) + 
                    Im(t(w_eight) %*% t(t(spec_mat) * sqrt(1 + Gamma * lambda))) %*% 
                    Im(t(t(t(spec_mat) * sqrt(1 + Gamma * lambda))) %*% t(des_mat)))
    scaler <- mean(diag(xtx))
    X_inv <- solve(xtx + scaler * reg_mat)
    bh <- as.vector(X_inv %*% (((t(Re(X_new) * weight_target)) %*% Gamma) - xtxy - scaler * reg_xtxy))
    b <- matrix(nrow = L, ncol = length(weight_h_exp[1, ]))
    bhh <- t(des_mat) %*% bh
    for (k in 1:L) {
        b[k, ] <- bhh[(k) + (0:(length(weight_h_exp[1, ]) - 1)) * L]
    }
    weight_cm <- matrix(w_eight, ncol = (length(weight_h_exp[1, ])))
    b <- b + weight_cm
    trth <- ((X_new) %*% (X_inv %*% t(Re(X_new)))) %*% (weight_target * Gamma)
    Proj_mat <- ((X_new) %*% (X_inv %*% t(Re(X_new))))
    res_mat <- diag(rep(1, dim(Proj_mat)[1])) - Proj_mat
    sum(abs(res_mat %*% (weight_target * Gamma))^2)
    resi <- res_mat %*% (weight_target * Gamma)
    t(Conj(resi)) %*% resi
    t((weight_target * Gamma)) %*% (t(Conj(res_mat)) %*% (res_mat)) %*% (weight_target * Gamma)
    degrees_freedom <- 2 * Re(sum(diag(t(Conj(res_mat)) %*% (res_mat)))) - 1
    Re(t(Conj(res_mat)) %*% (res_mat)) - Re(res_mat)
    freezed_degrees <- 2 * K + 1 - degrees_freedom
    2 * Re(sum(diag(Proj_mat)))
    sum(abs(Gamma * weight_target - trth)^2)
    trffkt <- matrix(nrow = K + 1, ncol = length(weight_h_exp[1, ]))
    trffkth <- trffkt
    trffkt[1, ] <- apply(b, 2, sum)
    trffkth[1, ] <- trffkt[1, ]
    for (j in 1:length(weight_h_exp[1, ])) {
        for (k in 0:(K)) {
            trffkt[k + 1, j] <- (b[, j] %*% exp((0+1i) * k * (0:(L - 1)) * pi/(K)))
        }
    }
    trt <- apply(((trffkt) * exp((0+1i) * (0 - Lag) * pi * (0:(K))/K)) * weight_h_exp, 1, sum)
    rever <- sum(abs(Gamma * weight_target - Re(trt) - (0+1i) * sqrt(1 + lambda * Gamma) * Im(trt))^2)/(2 * (K + 1)^2)
    MS_error <- sum((abs(Gamma * weight_target - trt)/expweight_vec)^2)/(2 * (K + 1)^2)
    Gamma_cp <- Gamma[1 + 0:as.integer(K * (cutoff/pi))]
    Gamma_cn <- Gamma[(2 + as.integer(K * (cutoff/pi))):(K + 1)]
    trt_cp <- (trt/expweight_vec)[1 + 0:as.integer(K * (cutoff/pi))]
    trt_cn <- (trt/expweight_vec)[(2 + as.integer(K * (cutoff/pi))):(K + 1)]
    weight_target_cp <- (weight_target/expweight_vec)[1 + 0:as.integer(K * (cutoff/pi))]
    weight_target_cn <- (weight_target/expweight_vec)[(2 + as.integer(K * (cutoff/pi))):(K + 1)]
    Accuracy <- sum(abs(Gamma_cp * weight_target_cp - abs(trt_cp))^2)/(2 * (K + 1)^2)
    Timeliness <- 4 * sum(abs(Gamma_cp) * abs(trt_cp) * sin(Arg(trt_cp)/2)^2 * weight_target_cp)/(2 * (K + 1)^2)
    Smoothness <- sum(abs(Gamma_cn * weight_target_cn - abs(trt_cn))^2)/(2 * (K + 1)^2)
    Shift_stopband <- 4 * sum(abs(Gamma_cn) * abs(trt_cn) * sin(Arg(trt_cn)/2)^2 * weight_target_cn)/(2 * (K + 1)^2)
    Accuracy + Timeliness + Smoothness + Shift_stopband - MS_error
    aic <- ifelse(degrees_freedom < K + 1 & degrees_freedom > 1, 
            log(rever) + 2 * (K - degrees_freedom + 1)/(degrees_freedom - 2), 
            NA)
    return(list(b = b, trffkt = trffkt, rever = rever, degrees_freedom = degrees_freedom, 
                    aic = aic, freezed_degrees = freezed_degrees, Accuracy = Accuracy, Smoothness = Smoothness, 
                    Timeliness = Timeliness, MS_error = MS_error))
}
MS_decomp_total <- function(Gamma, trffkt, weight_func, cutoff, Lag) {
    if (!(length(trffkt[, 1]) == length(weight_func[, 1]))) {
        len_w <- min(length(trffkt[, 1]), length(weight_func[, 1]))
        if (length(trffkt[, 1]) < length(weight_func[, 1])) {
            len_r <- (length(weight_func[, 1]) - 1)/(length(trffkt[, 1]) - 1)
            weight_funch <- weight_func[c(1, (1:(len_w - 1)) * len_r), ]
            trffkth <- trffkt
        }
        else {
            len_r <- 1/((length(weight_func[, 1]) - 1)/(length(trffkt[, 1]) - 1))
            trffkth <- trffkt[c(1, (1:(len_w - 1)) * len_r), ]
            weight_funch <- weight_func
        }
    }
    else {
        len_w <- length(trffkt[, 1])
        weight_funch <- weight_func
        trffkth <- trffkt
        Gammah <- Gamma
    }
    if (length(Gamma) > len_w) {
        len_r <- (length(Gamma) - 1)/(len_w - 1)
        Gammah <- Gamma[c(1, (1:(len_w - 1)) * len_r)]
    }
    weight_h <- weight_funch
    K <- length(weight_funch[, 1]) - 1
    weight_target <- weight_h[, 1]
    weight_h <- weight_h * exp(-(0+1i) * Arg(weight_target))
    weight_target <- Re(weight_target * exp(-(0+1i) * Arg(weight_target)))
    weight_h_exp <- as.matrix(weight_h[, 2:(dim(weight_h)[2])])
    trt <- apply(((trffkth) * exp((0+1i) * (0 - Lag) * pi * (0:(K))/K)) * weight_h_exp, 1, sum)
    MS_error <- sum((abs(Gammah * weight_target - trt))^2)/(2 * (K + 1)^2)
    Gamma_cp <- Gammah[1 + 0:as.integer(K * (cutoff/pi))]
    Gamma_cn <- Gammah[(2 + as.integer(K * (cutoff/pi))):(K + 1)]
    trt_cp <- trt[1 + 0:as.integer(K * (cutoff/pi))]
    trt_cn <- trt[(2 + as.integer(K * (cutoff/pi))):(K + 1)]
    weight_target_cp <- weight_target[1 + 0:as.integer(K * (cutoff/pi))]
    weight_target_cn <- weight_target[(2 + as.integer(K * (cutoff/pi))):(K + 1)]
    Accuracy <- sum(abs(Gamma_cp * weight_target_cp - abs(trt_cp))^2)/(2 * (K + 1)^2)
    Timeliness <- 4 * sum(abs(Gamma_cp) * abs(trt_cp) * sin(Arg(trt_cp)/2)^2 * weight_target_cp)/(2 * (K + 1)^2)
    Smoothness <- sum(abs(Gamma_cn * weight_target_cn - abs(trt_cn))^2)/(2 * (K + 1)^2)
    Shift_stopband <- 4 * sum(abs(Gamma_cn) * abs(trt_cn) * sin(Arg(trt_cn)/2)^2 *  weight_target_cn)/(2 * (K + 1)^2)
    Accuracy + Timeliness + Smoothness + Shift_stopband - MS_error
    return(list(Accuracy = Accuracy, Smoothness = Smoothness, Timeliness = Timeliness, MS_error = MS_error))
}
