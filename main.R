library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(broom)
library(dplyr)

plot_estimates <- function(fit, stan_data) {
    tidy(fit, pars = "theta", conf.int = TRUE) %>%
        bind_cols(naive = tapply(
            stan_data$y, stan_data$grp, function(x) mean(log(x)))) %>%
        ggplot(aes(estimate, term)) +
        geom_linerange(
            aes(xmin = conf.low, xmax = conf.high), colour = "orange", size = 3) +
        geom_point(shape = 21, size = 2.5, fill = "white") +
        geom_point(aes(naive, term), shape = 21, size = 4) +
        geom_vline(
            xintercept = mean(tapply(
                stan_data$y, stan_data$grp, function(x) mean(log(x)))),
            linetype = 3)
}


# Data
N <- c(5, 10, 20)
mu <- c(1, 5, 2.7)
set.seed(2020)
grp <- rep(seq_along(N), times = N)
y <- unlist(Map(function(n, meanlog) rlnorm(n, meanlog), N, mu))
stan_data <- list(N = sum(N), J = length(N), y = y, grp = grp)


# 1. Centred parametrisation
mod1 <- stan_model("cp_lognormal.stan")
fit1 <- sampling(object = mod1, data = stan_data, iter = 4000, seed = 2020)
summary(fit1, pars = c("theta", "mu_theta", "sigma_theta"))$summary
#                 mean     se_mean        sd       2.5%         25%       50%
#theta[1]    0.2824138 0.006508752 0.5826241 -0.8473542 -0.09846581 0.2754291
#theta[2]    5.3376859 0.004699512 0.4041809  4.5187907  5.07297626 5.3413570
#theta[3]    3.0580520 0.002967625 0.2864168  2.4891563  2.86236885 3.0576336
#mu_theta    2.9246911 0.026373804 1.8648984 -0.9152184  1.86809094 2.9457621
#sigma_theta 3.3292770 0.032916514 2.1642098  1.2496913  2.05277046 2.7801046
#                  75%    97.5%    n_eff      Rhat
#theta[1]    0.6616617 1.451716 8012.749 1.0000274
#theta[2]    5.6124935 6.102726 7396.835 0.9998538
#theta[3]    3.2534492 3.616276 9314.916 1.0000188
#mu_theta    3.9611710 6.741464 4999.940 1.0005051
#sigma_theta 3.9372162 8.707484 4322.859 0.9999944
plot_estimates(fit1, stan_data)
png("pairs_fit1.png", height = 6, width = 7, units = "in", res = 150)
pairs(fit1, pars = c("theta", "mu_theta", "log_sigma_theta"))
dev.off()


# 2. Non-centred parametrisation with Cauchy priors
mod2 <- stan_model("ncp_lognormal.stan")
fit2 <- sampling(object = mod2, data = stan_data, iter = 4000, seed = 2020)
summary(fit2, pars = c("theta", "mu_theta", "sigma_theta"))$summary
#                 mean     se_mean        sd       2.5%         25%       50%
#theta[1]    0.2874341 0.006494988 0.5736216 -0.8284460 -0.08609789 0.2759411
#theta[2]    5.3220404 0.004391849 0.4093376  4.4977139  5.04868842 5.3249077
#theta[3]    3.0644969 0.003102806 0.2854228  2.5044668  2.87820501 3.0633092
#mu_theta    3.0200486 0.039769094 1.7123739 -0.4843552  1.99999219 3.0028987
#sigma_theta 3.1652291 0.040395884 1.6448343  1.2707467  2.03742385 2.7391128
#                  75%    97.5%    n_eff     Rhat
#theta[1]    0.6593399 1.421572 7799.993 1.000028
#theta[2]    5.5914960 6.121670 8686.974 1.000006
#theta[3]    3.2534330 3.622223 8461.901 1.000132
#mu_theta    4.0012489 6.581174 1853.983 1.000884
#sigma_theta 3.8394788 7.619906 1657.945 1.007571
plot_estimates(fit2, stan_data)
png("pairs_fit2.png", height = 6, width = 7, units = "in", res = 150)
pairs(
    fit2,
    pars = c("theta", "mu_theta", "log_sigma_theta"),
    main = "Non-centred parametrisation, Cauchy priors on sigma's"))
dev.off()
png("pairs_raw_fit2.png", height = 6, width = 7, units = "in", res = 150)
pairs(
    fit2,
    pars = c("theta_raw", "sigma", "mu_theta", "log_sigma_theta"),
    main = "Non-centred parametrisation, Cauchy priors on sigma's")
dev.off()


# 3. Non-centred parametrisation with normal priors
mod3 <- stan_model("ncp_lognormal_normalpriors.stan")
fit3 <- sampling(object = mod3, data = stan_data, iter = 4000, seed = 2020)
summary(fit3, pars = c("theta", "mu_theta", "sigma_theta"))$summary
#                 mean     se_mean        sd       2.5%        25%       50%
#theta[1]    0.2582539 0.006465558 0.5822940 -0.8632974 -0.1304815 0.2564498
#theta[2]    5.3509043 0.004298230 0.4017381  4.5398392  5.0859275 5.3581667
#theta[3]    3.0689332 0.003085329 0.2865536  2.5025138  2.8811345 3.0721496
#mu_theta    2.9579136 0.046782353 1.9622719 -1.2196090  1.8440616 2.9489767
#sigma_theta 3.5846586 0.039733698 1.7723845  1.3619671  2.3031913 3.1730072
#                  75%    97.5%    n_eff      Rhat
#theta[1]    0.6380184 1.436320 8110.966 0.9996775
#theta[2]    5.6163492 6.134246 8735.879 0.9998688
#theta[3]    3.2603468 3.633102 8625.987 1.0001284
#mu_theta    4.0679984 7.064923 1759.358 1.0016995
#sigma_theta 4.4336332 8.214496 1989.747 1.0023587
plot_estimates(fit3, stan_data)
png("pairs_fit3.png", height = 6, width = 7, units = "in", res = 150)
pairs(
    fit3,
    pars = c("theta", "mu_theta", "log_sigma_theta"),
    main = "Non-centred parametrisation, normal priors on sigma's")
dev.off()
png("pairs_raw_fit3.png", height = 6, width = 7, units = "in", res = 150)
pairs(
    fit3,
    pars = c("theta_raw", "sigma", "mu_theta", "log_sigma_theta"),
    main = "Non-centred parametrisation, normal priors on sigma's")
dev.off()


# 4. Model in `brms`
library(brms)
fit3 <- brm(
    y ~  1 | grp,
    family = lognormal(),
    data = data.frame(y = y, grp = grp),
    seed = 2020)
fixef(fit3)[, "Estimate"] + ranef(fit3)$grp[, "Estimate", 1]
#        1         2         3
#0.2416648 5.3570510 3.0643358
