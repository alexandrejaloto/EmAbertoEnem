library (data.table)
library(tidyverse)
library (mirt)

banco = fread(
  'D:/Microdados/2018/MICRODADOS_ENEM_2018.csv',
  select = c(
    'NU_INSCRICAO',
    'Q001',
    'Q002',
    'Q003',
    'Q004',
    'Q006'
  )
)

# quem não sabe a escolaridade dos pais é NA
banco$Q001 = na_if(banco$Q001, 'H')
banco$Q002 = na_if(banco$Q002, 'H')

# quem não sabe a ocupação dos pais é NA
banco$Q003 = na_if(banco$Q003, 'F')
banco$Q004 = na_if(banco$Q004, 'F')

# recodificar a escolaridade para número
# descricao = LETTERS
nomes = vector ('list', length(LETTERS))
nomes[1:length(nomes)] = c(1:length(nomes))
names (nomes) = LETTERS
banco$Q001 = dplyr::recode(banco$Q001, !!!nomes)
banco$Q002 = dplyr::recode(banco$Q002, !!!nomes)
# ocupação
banco$Q003 = dplyr::recode(banco$Q003, !!!nomes)
banco$Q004 = dplyr::recode(banco$Q004, !!!nomes)
# e a renda (vai ser ordinal)
banco$Q006 = dplyr::recode(banco$Q006, !!!nomes)

calib = mirt (banco[,2:6], 1, TOL = .01, itemtype = 'graded')

nse = fscores(calib, response.pattern = banco[,2:6])
nse = data.frame (nse)

banco$nse = nse$F1

fwrite(banco, 'data-raw/nse_5itens.csv', sep = ';', dec = ',', row.names = FALSE, col.names = TRUE)
