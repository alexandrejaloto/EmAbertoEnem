library (data.table)
library (dplyr)
library (lmerTest)
library (stringr)

# rode os comandos do arquivo nse.R antes dos comandos abaixo

# caso já tenha rodado e gerado o arquivo nse_5itens.csv na pasta data-raw,
# passe para os comandos a seguir

source('R/funcao_analise_nse_aluno.R')

# as análises serão salvas na pasta data-raw
# para rodar linha a linha, abra o arquivo funcao_analise_nse_aluno.R

analise.nse.aluno('CH')

analise.nse.aluno('LC')

analise.nse.aluno('CN')

analise.nse.aluno('MT')

