source('Simulations_code.R')

seed.temp = 15

########################### von Mises noise ###########################

########################### S3 ###########################

S3 <- rep(0, 200)
S3_sim_FDR05 <- sim.study(S3, true.cpt = c(0), kappa = 2, seed = seed.temp, m = 100, FDR = 0.05)
S3_sim_FDR01 <- sim.study(S3, true.cpt = c(0), kappa = 2, seed = seed.temp, m = 100, FDR = 0.01)

########################### S4 ###########################

S4 <- c(rep(0, 50), rep(pi, 50))
S4_sim_kappa8_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 8, seed = seed.temp, m = 100, FDR = 0.05)
S4_sim_kappa8_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 8, seed = seed.temp, m = 100, FDR = 0.01)
S4_sim_kappa4_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 4, seed = seed.temp, m = 100, FDR = 0.05)
S4_sim_kappa4_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 4, seed = seed.temp, m = 100, FDR = 0.01)
S4_sim_kappa2_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 2, seed = seed.temp, m = 100, FDR = 0.05)
S4_sim_kappa2_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 2, seed = seed.temp, m = 100, FDR = 0.01)
S4_sim_kappa1_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 1, seed = seed.temp, m = 100, FDR = 0.05)
S4_sim_kappa1_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 1, seed = seed.temp, m = 100, FDR = 0.01)

########################### S5 ###########################

S5 <- c(rep(0, 50), rep(pi, 50), rep(1, 100))
S5_sim_kappa8_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 8, seed = seed.temp, m = 100, FDR = 0.05)
S5_sim_kappa8_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 8, seed = seed.temp, m = 100, FDR = 0.01)
S5_sim_kappa4_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 4, seed = seed.temp, m = 100, FDR = 0.05)
S5_sim_kappa4_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 4, seed = seed.temp, m = 100, FDR = 0.01)
S5_sim_kappa2_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 2, seed = seed.temp, m = 100, FDR = 0.05)
S5_sim_kappa2_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 2, seed = seed.temp, m = 100, FDR = 0.01)
S5_sim_kappa1_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 1, seed = seed.temp, m = 100, FDR = 0.05)
S5_sim_kappa1_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 1, seed = seed.temp, m = 100, FDR = 0.01)

S5_sim_kappa8_FDR01_B10000 <- sim.study(S5, true.cpt = c(50, 100), kappa = 8, seed = seed.temp, m = 100, FDR = 0.01, override_default = T)
S5_sim_kappa4_FDR01_B10000 <- sim.study(S5, true.cpt = c(50, 100), kappa = 4, seed = seed.temp, m = 100, FDR = 0.01, override_default = T)
S5_sim_kappa2_FDR01_B10000 <- sim.study(S5, true.cpt = c(50, 100), kappa = 2, seed = seed.temp, m = 100, FDR = 0.01, override_default = T)
S5_sim_kappa1_FDR01_B10000 <- sim.study(S5, true.cpt = c(50, 100), kappa = 1, seed = seed.temp, m = 100, FDR = 0.01, override_default = T)


########################### S6 ###########################

S6 <- c(rep(0, 30), rep(2, 30), rep(4, 30), rep(6, 30), rep(4, 30), rep(2, 30), rep(0, 30))
S6_sim_kappa8_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 8, seed = seed.temp, m = 100, FDR = 0.05)
S6_sim_kappa8_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 8, seed = seed.temp, m = 100, FDR = 0.01)
S6_sim_kappa4_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 4, seed = seed.temp, m = 100, FDR = 0.05)
S6_sim_kappa4_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 4, seed = seed.temp, m = 100, FDR = 0.01)
S6_sim_kappa2_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 2, seed = seed.temp, m = 100, FDR = 0.05)
S6_sim_kappa2_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 2, seed = seed.temp, m = 100, FDR = 0.01)
S6_sim_kappa1_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 1, seed = seed.temp, m = 100, FDR = 0.05)
S6_sim_kappa1_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 1, seed = seed.temp, m = 100, FDR = 0.01)

S6_sim_kappa8_FDR01_B10000 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 8, seed = seed.temp, m = 100, FDR = 0.01, override_default = T)
S6_sim_kappa4_FDR01_B10000 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 4, seed = seed.temp, m = 100, FDR = 0.01, override_default = T)
S6_sim_kappa2_FDR01_B10000 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 2, seed = seed.temp, m = 100, FDR = 0.01, override_default = T)
S6_sim_kappa1_FDR01_B10000 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 1, seed = seed.temp, m = 100, FDR = 0.01, override_default = T)


########################### S7 ###########################

S7 <- c(rep(1.5, 60), rep(3.3, 40), rep(5.2, 30), rep(1.5, 20))
S7_sim_kappa8_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 8, seed = seed.temp, m = 100, FDR = 0.05)
S7_sim_kappa8_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 8, seed = seed.temp, m = 100, FDR = 0.01)
S7_sim_kappa4_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 4, seed = seed.temp, m = 100, FDR = 0.05)
S7_sim_kappa4_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 4, seed = seed.temp, m = 100, FDR = 0.01)
S7_sim_kappa2_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 2, seed = seed.temp, m = 100, FDR = 0.05)
S7_sim_kappa2_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 2, seed = seed.temp, m = 100, FDR = 0.01)
S7_sim_kappa1_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 1, seed = seed.temp, m = 100, FDR = 0.05)
S7_sim_kappa1_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 1, seed = seed.temp, m = 100, FDR = 0.01)


########################### S8 ###########################

S8 <- c(rep(1, 150), rep(4, 150), rep(2, 200), rep(5, 100))
S8_sim_kappa8_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 8, seed = seed.temp, m = 100, FDR = 0.05)
S8_sim_kappa8_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 8, seed = seed.temp, m = 100, FDR = 0.01)
S8_sim_kappa4_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 4, seed = seed.temp, m = 100, FDR = 0.05)
S8_sim_kappa4_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 4, seed = seed.temp, m = 100, FDR = 0.01)
S8_sim_kappa2_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 2, seed = seed.temp, m = 100, FDR = 0.05)
S8_sim_kappa2_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 2, seed = seed.temp, m = 100, FDR = 0.01)
S8_sim_kappa1_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 1, seed = seed.temp, m = 100, FDR = 0.05)
S8_sim_kappa1_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 1, seed = seed.temp, m = 100, FDR = 0.01)


########################### wrapped Cauchy noise ###########################

########################### S3 ###########################

S3_sim_wrpCauchy_kappa0_86_FDR05 <- sim.study(S3, true.cpt = c(0), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S3_sim_wrpCauchy_kappa0_86_FDR01 <- sim.study(S3, true.cpt = c(0), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)


########################### S4 ###########################

S4_sim__wrpCauchy_kappa0_94_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S4_sim__wrpCauchy_kappa0_94_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S4_sim__wrpCauchy_kappa0_86_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S4_sim__wrpCauchy_kappa0_86_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S4_sim__wrpCauchy_kappa0_7_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S4_sim__wrpCauchy_kappa0_7_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S4_sim__wrpCauchy_kappa0_45_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S4_sim__wrpCauchy_kappa0_45_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)


########################### S5 ###########################

S5_sim__wrpCauchy_kappa0_94_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S5_sim__wrpCauchy_kappa0_94_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S5_sim__wrpCauchy_kappa0_86_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S5_sim__wrpCauchy_kappa0_86_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S5_sim__wrpCauchy_kappa0_7_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S5_sim__wrpCauchy_kappa0_7_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S5_sim__wrpCauchy_kappa0_45_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S5_sim__wrpCauchy_kappa0_45_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)


########################### S6 ###########################

S6_sim__wrpCauchy_kappa0_94_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S6_sim__wrpCauchy_kappa0_94_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S6_sim__wrpCauchy_kappa0_86_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S6_sim__wrpCauchy_kappa0_86_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S6_sim__wrpCauchy_kappa0_7_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S6_sim__wrpCauchy_kappa0_7_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S6_sim__wrpCauchy_kappa0_45_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S6_sim__wrpCauchy_kappa0_45_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)


########################### S7 ###########################

S7_sim__wrpCauchy_kappa0_94_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S7_sim__wrpCauchy_kappa0_94_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S7_sim__wrpCauchy_kappa0_86_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S7_sim__wrpCauchy_kappa0_86_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S7_sim__wrpCauchy_kappa0_7_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S7_sim__wrpCauchy_kappa0_7_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S7_sim__wrpCauchy_kappa0_45_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S7_sim__wrpCauchy_kappa0_45_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)


########################### S8 ###########################

S8_sim__wrpCauchy_kappa0_94_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S8_sim__wrpCauchy_kappa0_94_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S8_sim__wrpCauchy_kappa0_86_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S8_sim__wrpCauchy_kappa0_86_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S8_sim__wrpCauchy_kappa0_7_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S8_sim__wrpCauchy_kappa0_7_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)
S8_sim__wrpCauchy_kappa0_45_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.05)
S8_sim__wrpCauchy_kappa0_45_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpCauchy', FDR = 0.01)


########################### wrapped Normal noise ###########################

########################### S3 ###########################

S3_sim_wrpNormal_kappa0_86_FDR05 <- sim.study(S3, true.cpt = c(0), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S3_sim_wrpNormal_kappa0_86_FDR01 <- sim.study(S3, true.cpt = c(0), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)


########################### S4 ###########################

S4_sim__wrpNormal_kappa0_94_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S4_sim__wrpNormal_kappa0_94_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S4_sim__wrpNormal_kappa0_86_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S4_sim__wrpNormal_kappa0_86_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S4_sim__wrpNormal_kappa0_7_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S4_sim__wrpNormal_kappa0_7_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S4_sim__wrpNormal_kappa0_45_FDR05 <- sim.study(S4, true.cpt = c(50), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S4_sim__wrpNormal_kappa0_45_FDR01 <- sim.study(S4, true.cpt = c(50), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)


########################### S5 ###########################

S5_sim__wrpNormal_kappa0_94_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S5_sim__wrpNormal_kappa0_94_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S5_sim__wrpNormal_kappa0_86_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S5_sim__wrpNormal_kappa0_86_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S5_sim__wrpNormal_kappa0_7_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S5_sim__wrpNormal_kappa0_7_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S5_sim__wrpNormal_kappa0_45_FDR05 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S5_sim__wrpNormal_kappa0_45_FDR01 <- sim.study(S5, true.cpt = c(50, 100), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)


########################### S6 ###########################

S6_sim__wrpNormal_kappa0_94_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S6_sim__wrpNormal_kappa0_94_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S6_sim__wrpNormal_kappa0_86_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S6_sim__wrpNormal_kappa0_86_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S6_sim__wrpNormal_kappa0_7_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S6_sim__wrpNormal_kappa0_7_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S6_sim__wrpNormal_kappa0_45_FDR05 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S6_sim__wrpNormal_kappa0_45_FDR01 <- sim.study(S6, true.cpt = seq(30, 180, by = 30), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)


########################### S7 ###########################

S7_sim__wrpNormal_kappa0_94_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S7_sim__wrpNormal_kappa0_94_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S7_sim__wrpNormal_kappa0_86_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S7_sim__wrpNormal_kappa0_86_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S7_sim__wrpNormal_kappa0_7_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S7_sim__wrpNormal_kappa0_7_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S7_sim__wrpNormal_kappa0_45_FDR05 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S7_sim__wrpNormal_kappa0_45_FDR01 <- sim.study(S7, true.cpt = c(60,100,130), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)


########################### S8 ###########################

S8_sim__wrpNormal_kappa0_94_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S8_sim__wrpNormal_kappa0_94_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.94, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S8_sim__wrpNormal_kappa0_86_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S8_sim__wrpNormal_kappa0_86_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.86, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S8_sim__wrpNormal_kappa0_7_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S8_sim__wrpNormal_kappa0_7_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.7, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)
S8_sim__wrpNormal_kappa0_45_FDR05 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.05)
S8_sim__wrpNormal_kappa0_45_FDR01 <- sim.study(S8, true.cpt = c(150,300,500), kappa = 0.45, seed = seed.temp, m = 100, noise = 'wrpNormal', FDR = 0.01)

