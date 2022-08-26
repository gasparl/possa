# library('POSSA')

test_that("get_p gives correct p for simple htests", {
    expect_equal(0.07939414, get_p(t.test(extra ~ group, data = sleep)), tolerance = 0.0001)
    expect_equal(0.005585477, get_p(prop.test(c(83, 90, 129, 70), c(86, 93, 136, 82))), tolerance = 0.0001)
})


# ANOVA

test_that("get_p gives correct p for ANOVAs (aov, aovlist)", {
    expect_equal(list(block=0.0159387902081939,N=0.00437181182579938,P=0.474904092674434,K=0.0287950535002326,"N:P"=0.263165282877168,"N:K"=0.168647878500493,"P:K"=0.862752085685407),
        get_p(aov(yield ~ block + N * P * K, npk)),
        tolerance = 0.0001
    )

    expect_equal(list(block=0.0129536460360352,N=0.00368367600025441,P=0.475636972047535,"N:P"=0.26284056316726,K=0.0271137998320003),
        get_p(aov(
            terms(yield ~ block + N * P + K, keep.order = TRUE), npk
        ))
    )
    expect_equal(list("Error:block"=list("N:P:K"=0.525236141197407),"Error:Within"=list(N=0.00437181182579937,P=0.474904092674435,K=0.0287950535002327,"N:P"=0.263165282877168,"N:K"=0.168647878500492,"P:K"=0.862752085685407)), get_p(aov(yield ~ N * P * K + Error(block), npk)))
})


# ANOVA sim and pow

sampleAOV = function(samp_size) {
  dat_h0_human = faux::rnorm_multi(
    n = samp_size,
    vars = 2,
    mu = 1,
    sd = 1.03,
    r = 0.8
  )
  dat_h0_robot = faux::rnorm_multi(
    n = samp_size,
    vars = 2,
    mu = 1,
    sd = 1.03,
    r = 0.8
  )
  dat_h1_human = faux::rnorm_multi(
    n = samp_size,
    vars = 2,
    mu = c(1.03, 1.41),
    sd = 1.03,
    r = 0.8
  )
  dat_h1_robot = faux::rnorm_multi(
    n = samp_size,
    vars = 2,
    mu = c(0.98, 1.01),
    sd = 1.03,
    r = 0.8
  )
  list(
    # human_cheerful
    grp_1_human_cheerful_h0 = dat_h0_human$X1,
    grp_1_human_cheerful_h1 = dat_h1_human$X1,
    # human_sad
    grp_1_human_sad_h0  = dat_h0_human$X2,
    grp_1_human_sad_h1 = dat_h1_human$X2,
    # robot_cheerful
    grp_2_robot_cheerful_h0 = dat_h0_robot$X1,
    grp_2_robot_cheerful_h1 = dat_h1_robot$X1,
    # robot_sad
    grp_2_robot_sad_h0 = dat_h0_robot$X2,
    grp_2_robot_sad_h1 = dat_h1_robot$X2
  )
}

testAOV = function(grp_1_human_cheerful_h0,
                   grp_1_human_cheerful_h1,
                   grp_1_human_sad_h0,
                   grp_1_human_sad_h1,
                   grp_2_robot_cheerful_h0,
                   grp_2_robot_cheerful_h1,
                   grp_2_robot_sad_h0,
                   grp_2_robot_sad_h1) {
  len_grp1 = length(grp_1_human_cheerful_h0)
  len_grp2 = length(grp_2_robot_cheerful_h0)
  raw_data = data.frame(
    obs = c(
      grp_1_human_cheerful_h0,
      grp_1_human_sad_h0,
      grp_2_robot_cheerful_h0,
      grp_2_robot_sad_h0
    ),
    id = c(paste0('g1_', c(
      1:len_grp1, 1:len_grp1
    )),
    paste0('g2_', c(
      1:len_grp2, 1:len_grp2
    ))),
    voice = c(rep('human', len_grp1 * 2), rep('robot', len_grp2 * 2)),
    emotion = c(
      rep('cheerful', len_grp1),
      rep('sad', len_grp1),
      rep('cheerful', len_grp2),
      rep('sad', len_grp2)
    )
  )
  raw_data <<-raw_data
  aov_h0 = get_p(aov(obs ~ voice * emotion + Error(id / emotion), data =
                         raw_data))
  raw_data$obs = c(
    grp_1_human_cheerful_h1,
    grp_1_human_sad_h1,
    grp_2_robot_cheerful_h1,
    grp_2_robot_sad_h1
  )
  aov_h1 = get_p(aov(obs ~ voice * emotion + Error(id / emotion), data =
                         raw_data))

  return(
    c(
      p_voice_h0 = aov_h0$`Error:id`$voice,
      p_voice_h1 = aov_h1$`Error:id`$voice,
      p_emo_h0 = aov_h0$`Error:id:emotion`$emotion,
      p_emo_h1 = aov_h1$`Error:id:emotion`$emotion,
      p_interact_h0 = aov_h0$`Error:id:emotion`$`voice:emotion`,
      p_interact_h1 = aov_h1$`Error:id:emotion`$`voice:emotion`,
      p_sad_rob_vs_hum_h0 = get_p(t.test(grp_1_human_sad_h0, grp_2_robot_sad_h0, var.equal = TRUE)),
      p_sad_rob_vs_hum_h1 = get_p(t.test(grp_1_human_sad_h1, grp_2_robot_sad_h1, var.equal = TRUE)),
      p_cheer_rob_vs_hum_h0 = get_p(t.test(
        grp_1_human_cheerful_h0,
        grp_2_robot_cheerful_h0,
        var.equal = TRUE
      )),
      p_cheer_rob_vs_hum_h1 = get_p(t.test(
        grp_1_human_cheerful_h1,
        grp_2_robot_cheerful_h1,
        var.equal = TRUE
      ))
    )
  )
}

# do.call(testAOV, sampleAOV(100))

df_ps_aov_seq = sim(
  fun_obs = sampleAOV,
  n_obs = c(20, 30, 40),
  fun_test = testAOV,
  n_iter = 15
)

pow_aov_seq = pow(df_ps_aov_seq, alpha_locals = NA)
pow_aov_seq_any = pow(df_ps_aov_seq, alpha_locals = NA, multi_logic_a = 'any')

test_that("simulated p values for AOV", {

  expect_equal(tail(df_ps_aov_seq), structure(list(.iter=c(14,14,14,15,15,15),.look=c(1,2,3,1,2,3),.n_total=c(40,60,80,40,60,80),grp_1=c(20,30,40,20,30,40),grp_2=c(20,30,40,20,30,40),p_voice_h0=c(0.921951511569286,0.434953224524603,0.620283941206736,0.811199351981229,0.920828499692444,0.739704597929236),p_voice_h1=c(0.219369922784531,0.206316423956336,0.275337802811994,0.899873271845673,0.879565180178776,0.605729490357638),p_emo_h0=c(0.697713171008282,0.443022276325301,0.433085466970059,0.665914450522135,0.642419643192394,0.369060992775061),p_emo_h1=c(0.0257876206184076,0.00587520797388434,0.00057485456552077,0.0693064667341359,0.0219112769158262,0.0135440978341541),p_interact_h0=c(0.442595165986856,0.890947032101987,0.805001307851508,0.134192464139287,0.460952783089838,0.291897134474251),p_interact_h1=c(0.0623056981145687,0.0203427853374419,0.00176845453967279,4.2235064414923e-05,7.57645688516833e-06,2.72137336766717e-06),p_sad_rob_vs_hum_h0=c(0.762534441723391,0.39037565942833,0.556351967370459,0.502037573255412,0.88822554707346,0.990975401495593),p_sad_rob_vs_hum_h1=c(0.072350275077053,0.0373528788407116,0.0260549925368666,0.240392997152922,0.113942536219778,0.357244991740298),p_cheer_rob_vs_hum_h0=c(0.936708539273815,0.506617188033247,0.703568700843058,0.770392568489653,0.716764594766118,0.48975940711515),p_cheer_rob_vs_hum_h1=c(0.574183088684196,0.68746174114211,0.989199001502228,0.185175424864017,0.211845506017095,0.0657342977322442)),row.names=c(42L,41L,40L,45L,44L,43L),class=c("possa_sim_df","data.frame")))

})


test_that("power output for varied factors", {

  expect_equal(pow_aov_seq$pow, structure(list(look=c("1","2","3","totals"),n=c(40,60,80,180),n_rate=c(0.5,0.75,1,2.25),iters_sign_h0=c(0,0,0,0),iters_sign_h1=c(0,0,15,15),iters_sign_p_voice_h0=c(0,0,0,0),iters_sign_p_voice_h1=c(0,0,1,1),iters_sign_p_emo_h0=c(0,0,0,0),iters_sign_p_emo_h1=c(0,0,13,13),iters_sign_p_interact_h0=c(0,0,0,0),iters_sign_p_interact_h1=c(0,0,8,8),iters_sign_p_sad_rob_vs_hum_h0=c(0,0,0,0),iters_sign_p_sad_rob_vs_hum_h1=c(0,0,4,4),iters_sign_p_cheer_rob_vs_hum_h0=c(0,0,0,0),iters_sign_p_cheer_rob_vs_hum_h1=c(0,0,0,0),iters_remain_h0=c(15,15,0,30),iters_remain_h1=c(15,15,0,30),iters_stopped_h0=c(0,0,15,15),iters_stopped_h1=c(0,0,15,15),alpha_local_p_voice=structure(c(0.0332291666666667,0.0332291666666667,0.0332291666666667,0.0996875),class=c("numeric","possa_p")),ratio_sign_p_voice_h0=c(0,0,0,0),ratio_sign_p_voice_h1=c(0,0,0.0666666666666667,0.0666666666666667),alpha_local_p_emo=structure(c(0.0332291666666667,0.0332291666666667,0.0332291666666667,0.0996875),class=c("numeric","possa_p")),ratio_sign_p_emo_h0=c(0,0,0,0),ratio_sign_p_emo_h1=c(0,0,0.866666666666667,0.866666666666667),alpha_local_p_interact=structure(c(0.0332291666666667,0.0332291666666667,0.0332291666666667,0.0996875),class=c("numeric","possa_p")),ratio_sign_p_interact_h0=c(0,0,0,0),ratio_sign_p_interact_h1=c(0,0,0.533333333333333,0.533333333333333),alpha_local_p_sad_rob_vs_hum=structure(c(0.0332291666666667,0.0332291666666667,0.0332291666666667,0.0996875),class=c("numeric","possa_p")),ratio_sign_p_sad_rob_vs_hum_h0=c(0,0,0,0),ratio_sign_p_sad_rob_vs_hum_h1=c(0,0,0.266666666666667,0.266666666666667),alpha_local_p_cheer_rob_vs_hum=structure(c(0.0332291666666667,0.0332291666666667,0.0332291666666667,0.0996875),class=c("numeric","possa_p")),ratio_sign_p_cheer_rob_vs_hum_h0=c(0,0,0,0),ratio_sign_p_cheer_rob_vs_hum_h1=c(0,0,0,0),ratio_stopped_h0=c(0,0,1,1),ratio_stopped_h1=c(0,0,1,1),ratio_combined_sign_h0=c(0,0,0,0),ratio_combined_sign_h1=c(0,0,1,1),ratio_remain_h0=c(1,1,0,2),ratio_remain_h1=c(1,1,0,2),n_avg_prop_0=c(0,0,80,80),n_avg_prop_1=c(0,0,80,80)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")))

  expect_equal(pow_aov_seq_any$pow, structure(list(look=c("1","2","3","totals"),n=c(40,60,80,180),n_rate=c(0.5,0.75,1,2.25),iters_sign_h0=c(0,0,0,0),iters_sign_h1=c(8,6,1,15),iters_sign_p_voice_h0=c(0,0,0,0),iters_sign_p_voice_h1=c(0,0,0,0),iters_sign_p_emo_h0=c(0,0,0,0),iters_sign_p_emo_h1=c(5,6,1,12),iters_sign_p_interact_h0=c(0,0,0,0),iters_sign_p_interact_h1=c(3,1,1,5),iters_sign_p_sad_rob_vs_hum_h0=c(0,0,0,0),iters_sign_p_sad_rob_vs_hum_h1=c(0,1,1,2),iters_sign_p_cheer_rob_vs_hum_h0=c(0,0,0,0),iters_sign_p_cheer_rob_vs_hum_h1=c(0,1,0,1),iters_remain_h0=c(15,15,0,30),iters_remain_h1=c(7,1,0,8),iters_stopped_h0=c(0,0,15,15),iters_stopped_h1=c(8,6,1,15),alpha_local_p_voice=structure(c(0.0223307291666667,0.0223307291666667,0.0223307291666667,0.0669921875),class=c("numeric","possa_p")),ratio_sign_p_voice_h0=c(0,0,0,0),ratio_sign_p_voice_h1=c(0,0,0,0),alpha_local_p_emo=structure(c(0.0223307291666667,0.0223307291666667,0.0223307291666667,0.0669921875),class=c("numeric","possa_p")),ratio_sign_p_emo_h0=c(0,0,0,0),ratio_sign_p_emo_h1=c(0.333333333333333,0.4,0.0666666666666667,0.8),alpha_local_p_interact=structure(c(0.0223307291666667,0.0223307291666667,0.0223307291666667,0.0669921875),class=c("numeric","possa_p")),ratio_sign_p_interact_h0=c(0,0,0,0),ratio_sign_p_interact_h1=c(0.2,0.0666666666666667,0.0666666666666667,0.333333333333333),alpha_local_p_sad_rob_vs_hum=structure(c(0.0223307291666667,0.0223307291666667,0.0223307291666667,0.0669921875),class=c("numeric","possa_p")),ratio_sign_p_sad_rob_vs_hum_h0=c(0,0,0,0),ratio_sign_p_sad_rob_vs_hum_h1=c(0,0.0666666666666667,0.0666666666666667,0.133333333333333),alpha_local_p_cheer_rob_vs_hum=structure(c(0.0223307291666667,0.0223307291666667,0.0223307291666667,0.0669921875),class=c("numeric","possa_p")),ratio_sign_p_cheer_rob_vs_hum_h0=c(0,0,0,0),ratio_sign_p_cheer_rob_vs_hum_h1=c(0,0.0666666666666667,0,0.0666666666666667),ratio_stopped_h0=c(0,0,1,1),ratio_stopped_h1=c(0.533333333333333,0.4,0.0666666666666667,1),ratio_combined_sign_h0=c(0,0,0,0),ratio_combined_sign_h1=c(0.533333333333333,0.4,0.0666666666666667,1),ratio_remain_h0=c(1,1,0,2),ratio_remain_h1=c(0.466666666666667,0.0666666666666667,0,0.533333333333333),n_avg_prop_0=c(0,0,80,80),n_avg_prop_1=c(21.3333333333333,24,5.33333333333333,50.6666666666667)),row.names=c(NA,-4L),class=c("possa_pow_df","data.frame")))

})

# cat('DF: ', gsub(" ", "", paste(deparse(dput(
#   pow_aov_seq_extra
# )), collapse = ''), fixed = TRUE))
#
