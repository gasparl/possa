library('POSSA')

# custom t-test function ####
ttest = function(x, y) {
  t_info = stats::t.test(x, y, paired = TRUE, var.equal = T)
  sdx = sd(x)
  sdy = sd(y)
  corr = cor(x, y)
  sd_p = sqrt((sdx ** 2 + sdy ** 2) - 2 * corr * sdx * sdy)
  return(list(
    pval = t_info$p.value,
    mean_diff = as.numeric(t_info$estimate),
    corr =  corr,
    smd = as.numeric(t_info$estimate) / sd_p
  ))
}

my_samp = function(samp_size, h1_mean, h1_corr) {
  correlated_samples = faux::rnorm_multi(
    n = samp_size,
    vars = 3,
    mu = c(0, 0, h1_mean),
    sd = 5,
    r = c(h1_corr, h1_corr, 0)
  )
  list(
    GRP_v1 = correlated_samples$X1, # correlated with both X2 and X3
    GRP_v2_h0 = correlated_samples$X2, # correlated only with X1
    GRP_v2_h1 = correlated_samples$X3  # correlated only with X1
  )
}

my_test = function(GRP_v1, GRP_v2_h0, GRP_v2_h1) {
  t0 = ttest(GRP_v2_h0, GRP_v1)
  t1 = ttest(GRP_v2_h1, GRP_v1)
  return(c(
    p_h0 = t0$pval,
    m_diff_0 = t0$mean_diff,
    corr_0 = t0$corr,
    smd_0 = t0$smd,
    p_h1 = t1$pval,
    m_diff_1 = t1$mean_diff,
    corr_1 = t1$corr,
    smd_1 = t1$smd
  ))
}

# do.call(my_test, my_samp(30, 1, .5))

# run simulation ####
# varied parameters
df_ps_facts = sim(
  fun_obs = list(
    my_samp,
    h1_mean = c(1.5, 2.5, 3.5),
    h1_corr = c(0, 0.5)
    #h1_corr = c(0.5)
  ),
  n_obs = c(30, 60, 90),
  fun_test = my_test,
  pair = TRUE, n_iter = 100
)

pow_base = pow(df_ps_facts)

pow_subset = pow(df_ps_facts[df_ps_facts$h1_mean == 1.5 & df_ps_facts$h1_corr == 0,])

pow_new_global = pow(df_ps_facts, alpha_global = .001)

pow_NA = pow(df_ps_facts, alpha_locals = NA)

pow_fut = pow(df_ps_facts, alpha_locals = NA, fut_locals = c(0.5, 0.5))

test_that("simulated p values for varied factors", {
  expect_equal(head(df_ps_facts), structure(list(.iter=c(1,1,1,2,2,2),.look=c(1,2,3,1,2,3),.n_total=c(30,60,90,30,60,90),h1_mean=c(1.5,1.5,1.5,1.5,1.5,1.5),h1_corr=c(0,0,0,0,0,0),GRP=c(30,60,90,30,60,90),p_h0=c(0.468608419576857,0.45897940029039,0.688133841266967,0.695013807147141,0.552719198901775,0.156907317765909),m_diff_0=c(0.90071285984012,0.680343605464961,-0.3017988843388,0.48658399713647,0.519071064987721,1.01397902691137),corr_0=c(-0.00203754829643567,0.0299124934340033,0.0906565855517084,0.162705951142637,0.0812203952370393,0.13097050577546),smd_0=c(0.134080530005479,0.0962328498243775,-0.04244841445377,0.072296826260271,0.0770864574123762,0.150481939968026),p_h1=c(0.2833345980455,0.232925171155621,0.133813023833551,0.25011908636024,0.088812209659603,0.00946901367179193),m_diff_1=c(1.35341392796442,1.21964326022608,1.26530374788781,1.72543120316463,1.68685613569459,2.01515052680961),corr_1=c(-0.164172352230674,-0.113146965180433,-0.0902595219283293,-0.211171197948254,-0.162336458358617,-0.042396180473315),smd_1=c(0.199577217289796,0.155594917868192,0.159488239249683,0.214262162442653,0.223375365974416,0.279562548114999)),row.names=c(3L,2L,1L,6L,5L,4L),class=c("possa_sim_df","data.frame")))


  expect_equal(tail(df_ps_facts), structure(list(.iter=c(99,99,99,100,100,100),.look=c(1,2,3,1,2,3),.n_total=c(30,60,90,30,60,90),h1_mean=c(3.5,3.5,3.5,3.5,3.5,3.5),h1_corr=c(0.5,0.5,0.5,0.5,0.5,0.5),GRP=c(30,60,90,30,60,90),p_h0=c(0.793299853051033,0.449389659955649,0.654748672998688,0.360168219372326,0.413635008290108,0.30222648939799),m_diff_0=c(0.239320724427256,0.468196375980935,0.229367208929758,-0.949757561885943,0.565265661999364,0.590693479369086),corr_0=c(0.44410559035086,0.567251090943552,0.574018223236826,0.241204523472943,0.468706935940721,0.371782140311933),smd_0=c(0.0482826463753128,0.0983096668771111,0.0472956661620324,-0.169751052147242,0.106291669171268,0.109381904622381),p_h1=c(0.000112025616318969,9.81743492951011e-10,5.97642974841143e-11,5.45973736971821e-06,1.64342239670229e-09,1.73367212848835e-09),m_diff_1=c(3.68519600144514,4.79505541192009,4.06963955350423,2.9152097751864,3.84357292039329,3.18271472517354),corr_1=c(0.59838689404913,0.561305566061908,0.55401893234951,0.732381578663456,0.621030649338328,0.544406732049736),smd_1=c(0.815085189632904,0.937249038667122,0.784436435865142,1.01387407896267,0.920201278966792,0.707464913054077)),row.names=c(1797L,1796L,1795L,1800L,1799L,1798L),class=c("possa_sim_df","data.frame")))
})


test_that("power output for varied factors", {
  expect_equal(pow_base[2:5],list(pow_1.5_0.5=structure(list(look=c("1","2","3","totals"),n=c(30,60,90,180),n_rate=c(0.333333333333333,0.666666666666667,1,2),iters_sign_h0=c(0,0,2,2),iters_sign_h1=c(0,0,83,83),iters_sign_p_h0=c(0,0,2,2),iters_sign_p_h1=c(0,0,83,83),iters_remain_h0=c(100,100,0,200),iters_remain_h1=c(100,100,0,200),iters_stopped_h0=c(0,0,100,100),iters_stopped_h1=c(0,0,100,100),alpha_local_p=structure(c(0,0,0.05,0.05),class=c("numeric","possa_p")),ratio_sign_p_h0=c(0,0,0.02,0.02),ratio_sign_p_h1=c(0,0,0.83,0.83),ratio_stopped_h0=c(0,0,1,1),ratio_stopped_h1=c(0,0,1,1),ratio_combined_sign_h0=c(0,0,0.02,0.02),ratio_combined_sign_h1=c(0,0,0.83,0.83),ratio_remain_h0=c(1,1,0,2),ratio_remain_h1=c(1,1,0,2),n_avg_prop_0=c(0,0,90,90),n_avg_prop_1=c(0,0,90,90)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")),pow_2.5_0=structure(list(look=c("1","2","3","totals"),n=c(30,60,90,180),n_rate=c(0.333333333333333,0.666666666666667,1,2),iters_sign_h0=c(0,0,7,7),iters_sign_h1=c(0,0,89,89),iters_sign_p_h0=c(0,0,7,7),iters_sign_p_h1=c(0,0,89,89),iters_remain_h0=c(100,100,0,200),iters_remain_h1=c(100,100,0,200),iters_stopped_h0=c(0,0,100,100),iters_stopped_h1=c(0,0,100,100),alpha_local_p=structure(c(0,0,0.05,0.05),class=c("numeric","possa_p")),ratio_sign_p_h0=c(0,0,0.07,0.07),ratio_sign_p_h1=c(0,0,0.89,0.89),ratio_stopped_h0=c(0,0,1,1),ratio_stopped_h1=c(0,0,1,1),ratio_combined_sign_h0=c(0,0,0.07,0.07),ratio_combined_sign_h1=c(0,0,0.89,0.89),ratio_remain_h0=c(1,1,0,2),ratio_remain_h1=c(1,1,0,2),n_avg_prop_0=c(0,0,90,90),n_avg_prop_1=c(0,0,90,90)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")),pow_2.5_0.5=structure(list(look=c("1","2","3","totals"),n=c(30,60,90,180),n_rate=c(0.333333333333333,0.666666666666667,1,2),iters_sign_h0=c(0,0,2,2),iters_sign_h1=c(0,0,100,100),iters_sign_p_h0=c(0,0,2,2),iters_sign_p_h1=c(0,0,100,100),iters_remain_h0=c(100,100,0,200),iters_remain_h1=c(100,100,0,200),iters_stopped_h0=c(0,0,100,100),iters_stopped_h1=c(0,0,100,100),alpha_local_p=structure(c(0,0,0.05,0.05),class=c("numeric","possa_p")),ratio_sign_p_h0=c(0,0,0.02,0.02),ratio_sign_p_h1=c(0,0,1,1),ratio_stopped_h0=c(0,0,1,1),ratio_stopped_h1=c(0,0,1,1),ratio_combined_sign_h0=c(0,0,0.02,0.02),ratio_combined_sign_h1=c(0,0,1,1),ratio_remain_h0=c(1,1,0,2),ratio_remain_h1=c(1,1,0,2),n_avg_prop_0=c(0,0,90,90),n_avg_prop_1=c(0,0,90,90)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")),pow_3.5_0=structure(list(look=c("1","2","3","totals"),n=c(30,60,90,180),n_rate=c(0.333333333333333,0.666666666666667,1,2),iters_sign_h0=c(0,0,8,8),iters_sign_h1=c(0,0,100,100),iters_sign_p_h0=c(0,0,8,8),iters_sign_p_h1=c(0,0,100,100),iters_remain_h0=c(100,100,0,200),iters_remain_h1=c(100,100,0,200),iters_stopped_h0=c(0,0,100,100),iters_stopped_h1=c(0,0,100,100),alpha_local_p=structure(c(0,0,0.05,0.05),class=c("numeric","possa_p")),ratio_sign_p_h0=c(0,0,0.08,0.08),ratio_sign_p_h1=c(0,0,1,1),ratio_stopped_h0=c(0,0,1,1),ratio_stopped_h1=c(0,0,1,1),ratio_combined_sign_h0=c(0,0,0.08,0.08),ratio_combined_sign_h1=c(0,0,1,1),ratio_remain_h0=c(1,1,0,2),ratio_remain_h1=c(1,1,0,2),n_avg_prop_0=c(0,0,90,90),n_avg_prop_1=c(0,0,90,90)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame"))))

  expect_equal(pow_subset$pow, structure(list(look=c("1","2","3","totals"),n=c(30,60,90,180),n_rate=c(0.333333333333333,0.666666666666667,1,2),iters_sign_h0=c(0,0,2,2),iters_sign_h1=c(0,0,51,51),iters_sign_p_h0=c(0,0,2,2),iters_sign_p_h1=c(0,0,51,51),iters_remain_h0=c(100,100,0,200),iters_remain_h1=c(100,100,0,200),iters_stopped_h0=c(0,0,100,100),iters_stopped_h1=c(0,0,100,100),alpha_local_p=structure(c(0,0,0.05,0.05),class=c("numeric","possa_p")),ratio_sign_p_h0=c(0,0,0.02,0.02),ratio_sign_p_h1=c(0,0,0.51,0.51),ratio_stopped_h0=c(0,0,1,1),ratio_stopped_h1=c(0,0,1,1),ratio_combined_sign_h0=c(0,0,0.02,0.02),ratio_combined_sign_h1=c(0,0,0.51,0.51),ratio_remain_h0=c(1,1,0,2),ratio_remain_h1=c(1,1,0,2),n_avg_prop_0=c(0,0,90,90),n_avg_prop_1=c(0,0,90,90)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")))

  expect_equal(pow_new_global$pow_2.5_0.5, structure(list(look=c("1","2","3","totals"),n=c(30,60,90,180),n_rate=c(0.333333333333333,0.666666666666667,1,2),iters_sign_h0=c(0,0,0,0),iters_sign_h1=c(0,0,97,97),iters_sign_p_h0=c(0,0,0,0),iters_sign_p_h1=c(0,0,97,97),iters_remain_h0=c(100,100,0,200),iters_remain_h1=c(100,100,0,200),iters_stopped_h0=c(0,0,100,100),iters_stopped_h1=c(0,0,100,100),alpha_local_p=structure(c(0,0,0.001,0.001),class=c("numeric","possa_p")),ratio_sign_p_h0=c(0,0,0,0),ratio_sign_p_h1=c(0,0,0.97,0.97),ratio_stopped_h0=c(0,0,1,1),ratio_stopped_h1=c(0,0,1,1),ratio_combined_sign_h0=c(0,0,0,0),ratio_combined_sign_h1=c(0,0,0.97,0.97),ratio_remain_h0=c(1,1,0,2),ratio_remain_h1=c(1,1,0,2),n_avg_prop_0=c(0,0,90,90),n_avg_prop_1=c(0,0,90,90)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")))

  expect_equal(pow_NA$pow_1.5_0.5, structure(list(look=c("1","2","3","totals"),n=c(30,60,90,180),n_rate=c(0.333333333333333,0.666666666666667,1,2),iters_sign_h0=c(2,2,1,5),iters_sign_h1=c(40,30,15,85),iters_sign_p_h0=c(2,2,1,5),iters_sign_p_h1=c(40,30,15,85),iters_remain_h0=c(98,96,0,194),iters_remain_h1=c(60,30,0,90),iters_stopped_h0=c(2,2,96,100),iters_stopped_h1=c(40,30,30,100),alpha_local_p=structure(c(0.0491666666666667,0.0491666666666667,0.0491666666666667,0.1475),class=c("numeric","possa_p")),ratio_sign_p_h0=c(0.02,0.02,0.01,0.05),ratio_sign_p_h1=c(0.4,0.3,0.15,0.85),ratio_stopped_h0=c(0.02,0.02,0.96,1),ratio_stopped_h1=c(0.4,0.3,0.3,1),ratio_combined_sign_h0=c(0.02,0.02,0.01,0.05),ratio_combined_sign_h1=c(0.4,0.3,0.15,0.85),ratio_remain_h0=c(0.98,0.96,0,1.94),ratio_remain_h1=c(0.6,0.3,0,0.9),n_avg_prop_0=c(0.6,1.2,86.4,88.2),n_avg_prop_1=c(12,18,27,57)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")))

  expect_equal(pow_fut$pow_3.5_0.5, structure(list(look=c("1","2","3","totals"),n=c(30,60,90,180),n_rate=c(0.333333333333333,0.666666666666667,1,2),iters_sign_h0=c(3,2,0,5),iters_sign_h1=c(98,2,0,100),iters_fut_h0=c(51,18,26,95),iters_fut_h1=c(0,0,0,0),iters_sign_p_h0=c(3,2,0,5),iters_sign_p_h1=c(98,2,0,100),iters_futil_p_h0=c(50,1,0,51),iters_futil_p_h1=c(0,0,0,0),iters_remain_h0=c(46,26,0,72),iters_remain_h1=c(2,0,0,2),iters_stopped_h0=c(54,20,26,100),iters_stopped_h1=c(98,2,0,100),alpha_local_p=structure(c(0.0266666666666667,0.0266666666666667,0.0266666666666667,0.08),class=c("numeric","possa_p")),futil_local_p=structure(c(0.5,0.5,NA,NA),class=c("numeric","possa_p","possa_futility")),ratio_sign_p_h0=c(0.03,0.02,0,0.05),ratio_futil_p_h0=c(0.5,0.01,0,0.51),ratio_sign_p_h1=c(0.98,0.02,0,1),ratio_futil_p_h1=c(0,0,0,0),ratio_stopped_h0=c(0.54,0.2,0.26,1),ratio_stopped_h1=c(0.98,0.02,0,1),ratio_combined_sign_h0=c(0.03,0.02,0,0.05),ratio_combined_sign_h1=c(0.98,0.02,0,1),ratio_combined_fut_h0=c(0.51,0.18,0.26,0.95),ratio_combined_fut_h1=c(0,0,0,0),ratio_remain_h0=c(0.46,0.26,0,0.72),ratio_remain_h1=c(0.02,0,0,0.02),n_avg_prop_0=c(16.2,12,23.4,51.6),n_avg_prop_1=c(29.4,1.2,0,30.6)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")))

})


# cat('DF: ', gsub(" ", "", paste(deparse(dput(
#   pow_fut$df_3.5_0.5
# )), collapse = ''), fixed = TRUE))

