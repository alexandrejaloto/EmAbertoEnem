analise.nse.aluno = function(disc)
{

  # função para padronizar os coeficientes (será usada mais adiante)
  lm.beta.lmer <- function(mod) {
    b <- fixef(mod)[-1]
    sd.x <- apply(getME(mod,"X")[,-1],2,sd)
    sd.y <- sd(getME(mod,"y"))
    b*sd.x/sd.y
  }

  # importar microdados
  bruto = fread (
    'D:/Microdados/2018/MICRODADOS_ENEM_2018.csv',
    select = c(
      paste0 (
        c('TP_PRESENCA_', 'NU_NOTA_'),
        rep (c(disc),each = 2)
      ),
      'NU_INSCRICAO',
      'SG_UF_RESIDENCIA',
      'NU_IDADE',
      'TP_SEXO',
      'TP_COR_RACA',
      'TP_ST_CONCLUSAO',
      'TP_ESCOLA',
      'TP_ENSINO',
      'CO_ESCOLA',
      'SG_UF_ESC',
      'TP_DEPENDENCIA_ADM_ESC',
      'TP_LOCALIZACAO_ESC',
      "CO_MUNICIPIO_ESC"
    )
  )

  # filtrar com os critérios de inclusão
  banco = dplyr::filter (
    bruto,
    # ensino regular ou EJA
    TP_ENSINO %in% c(1,3) &
      # 2 == Estou cursando e concluirei o Ensino Médio em 2018
      TP_ST_CONCLUSAO == 2 &
      NU_IDADE >= 16 &
      get (paste0('TP_PRESENCA_', disc)) == 1
  )

  # importar dados do NSE
  nse = fread ('data-raw/nse_5itens.csv', sep = ';', dec = ',')

  # juntar as bases
  banco = left_join(banco, nse[,c('NU_INSCRICAO', 'nse')], by = 'NU_INSCRICAO')

  # retirar quem é NA (mesmo que não tire aqui, na função ele vai ser excluído. então já tira agora para ficar mais rápido)
  banco = banco %>%
    dplyr::filter (
      !is.na(TP_SEXO) &
        !is.na(TP_COR_RACA) &
        # quem não declarou cor também é NA
        TP_COR_RACA != 0 &
        !is.na(nse) &
        !is.na(TP_DEPENDENCIA_ADM_ESC)
    )

  # como uma mesma escola pode ter diferentes modalidades
  # de ensino, vamos fazer como se cada modalidade fosse um
  # código de escola diferente
  banco$cod.escola = stringr::str_c(banco$CO_ESCOLA, banco$TP_ENSINO)

  # qual é a frequência de escola?
  freq = table (banco$cod.escola)

  # que escolas têm o mínimo?
  min = names (which (freq >= 30))

  # selecionar os alunos
  banco = filter(banco, banco$cod.escola %in% min)

  # quem tem código de escola
  banco = banco [!is.na(banco$CO_ESCOLA),]

  # alterar para fator as variáveis
  banco$TP_ST_CONCLUSAO = as.factor(banco$TP_ST_CONCLUSAO)
  banco$TP_ENSINO = as.factor(banco$TP_ENSINO)
  banco$TP_DEPENDENCIA_ADM_ESC = as.factor(banco$TP_DEPENDENCIA_ADM_ESC)
  banco$TP_COR_RACA = as.factor(banco$TP_COR_RACA)

  # variável atraso escolar
  banco$ATRASO = ifelse(banco$NU_IDADE > 18, 1, 0)

  # bora deixar as variáveis em dummy
  banco$BRANCA = ifelse(banco$TP_COR_RACA == 1, 1, 0)
  banco$FEDERAL = ifelse(banco$TP_DEPENDENCIA_ADM_ESC == 1, 1, 0)
  banco$ESTADUAL = ifelse(banco$TP_DEPENDENCIA_ADM_ESC == 2, 1, 0)
  banco$PRIVADA = ifelse(banco$TP_DEPENDENCIA_ADM_ESC == 4, 1, 0)


  # modelo nulo
  form = as.formula(paste0(
    paste0('NU_NOTA_',disc),
    '~ 1 + (1|cod.escola)'))
  mod0 = lmer (
    form,
    REML = FALSE,
    control=lmerControl(optimizer="bobyqa",
                        optCtrl=list(maxfun=2e5)),
    data=banco
  )

  sum0 = summary(mod0)

  # caso não convirja, insere valores iniciais para facilitar convergência
  i0 = 0
  while (length (sum0$optinfo$conv$lme4$messages) > 0)
  {
    i0 = i0+1
    ss = getME(mod0 ,c("theta","fixef"))
    mod0 = update(
      mod0, start=ss,
      control=lmerControl(optCtrl=list(maxfun=2e4)))
    sum0 = summary(mod0)
  }

  # salvar coeficientes
  sum0$coefficients %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_coef_0.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar AIC
  sum0$AICtab %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_ajuste_0.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar variâncias e covariâncias
  sum0$varcor %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_var_0.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # variáveis do 1o nível
  form = as.formula(
    paste0(
      paste0('NU_NOTA_',disc),
      '~ 1 +
  ATRASO +
  TP_SEXO +
  BRANCA +
  nse +
  (1|cod.escola)'
    )
  )

  # modelo com variáveis de 1o nível
  mod1 = lmer (
    form,
    REML = FALSE,
    control=lmerControl(optimizer="bobyqa",
                        optCtrl=list(maxfun=2e5)),
    data=banco
  )

  sum1 = summary(mod1)

  # caso não convirja
  i1 = 0
  while (length (sum1$optinfo$conv$lme4$messages) > 0)
  {
    i1 = i1+1
    ss = getME(mod1 ,c("theta","fixef"))
    mod1 = update(
      mod1, start=ss,
      control=lmerControl(optCtrl=list(maxfun=2e4)))
    sum1 = summary(mod1)
  }

  # salvar coeficientes
  sum1$coefficients %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_coef_1.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar AIC
  sum1$AICtab %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_ajuste_1.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar variâncias e covariâncias
  sum1$varcor %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_var_1.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar coeficientes padronizados
  lm.beta.lmer(mod1) %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_beta_1.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # variáveis do 2o nível
  form = as.formula(
    paste0(
      paste0('NU_NOTA_',disc),
      '~ 1 +
  ATRASO +
  TP_SEXO +
  BRANCA +
  nse +
  ESTADUAL + PRIVADA + FEDERAL +
  (1|cod.escola)'
    )
  )

  # modelo com variáveis do 2o nível
  mod2.1 = lmer (
    form,
    REML = FALSE,
    control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
    data=banco
  )

  sum2.1 = summary(mod2.1)

  # caso não convirja
  i2.1 = 0
  while (length (sum2.1$optinfo$conv$lme4$messages) > 0)
  {
    i2.1 = i2.1+1
    ss = getME(mod2.1 ,c("theta","fixef"))
    mod2.1 = update(
      mod2.1, start=ss,
      control=lmerControl(optCtrl=list(maxfun=2e4)))
    sum2.1 = summary(mod2.1)
  }

  # salvar coeficientes
  sum2.1$coefficients %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_coef_2_1.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar AIC
  sum2.1$AICtab %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_ajuste_2_1.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar variâncias e covariâncias
  sum2.1$varcor %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_var_2_1.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # padronizar coeficientes
  lm.beta.lmer(mod2.1) %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_beta_2_1.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # variáveis do 2o nível
  # nse com slope randômico
  form = as.formula(
    paste0(
      paste0('NU_NOTA_',disc),
      '~ 1 +
  ATRASO +
  TP_SEXO +
  BRANCA +
  nse +
  ESTADUAL + PRIVADA + FEDERAL +
  (nse|cod.escola)'
    )
  )

  # modelo com slope randômico para NSE
  mod2.2 = lmer (
    form,
    REML = FALSE,
    control=lmerControl(optimizer="bobyqa",
                        optCtrl=list(maxfun=2e5)),
    data=banco
  )

  sum2.2 = summary (mod2.2)

  # caso não convirja
  i2.2 = 0
  while (length (sum2.2$optinfo$conv$lme4$messages) > 0)
  {
    i2.2 = i2.2+1
    ss = getME(mod2.2 ,c("theta","fixef"))
    mod2.2 = update(
      mod2.2, start=ss,
    )
    sum2.2 = summary(mod2.2)
  }

  # salvar coeficientes
  sum2.2$coefficients %>%
    data.frame() %>%
    fwrite(paste0 ('data-raw/', disc, '_coef_2_2.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar AIC
  sum2.2$AICtab %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_ajuste_2_2.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar variâncias e covariâncias
  sum2.2$varcor %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_var_2_2.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # padronizar coeficientes
  lm.beta.lmer(mod2.2) %>%
    data.frame() %>%
    fwrite(paste0('data-raw/', disc, '_beta_2_2.csv'),
           row.names = TRUE,
           sep = ';', dec = ',')

  # salvar informações sobre a convergência. se salvar o valor 0, significa que convergiu de primeira.
  # se salvar com valor 1, convergiu na 2a iteração.
  cat(paste(
    paste0('mod0 = ', i0),
    paste0('mod1 = ', i1),
    paste0('mod2.1 = ', i2.1),
    paste0('mod2.2 = ', i2.2),
    sep = '\n'),
    file = paste0('data-raw/', disc, '_convergencia.txt'))

  # salvar o banco utilizado
  fwrite(
    banco,
    paste0('data-raw/', disc, '_banco.csv'),
    sep = ';',
    dec = ',')

  # salvar RData com todas as informações
  save(
    mod0, mod1, mod2.1, mod2.2,
    file = paste0('data-raw/', disc, '_mod.RData'))
}
