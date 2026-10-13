# =============================================================================
# RELATORIO DE REVISAO DOS SCRIPTS  -  arbovirus-surveillance-brazil-r
# -----------------------------------------------------------------------------
# Autora da revisao: Claude (Cowork)  |  Solicitado por: Sofia Lirio (LII LAB)
# Data: 2026-10-13
#
# Este arquivo e um RELATORIO (so comentarios - abre no Positron sem executar
# nada). Ele lista, script por script:
#   (A) as ALTERACOES que foram efetivamente aplicadas nos arquivos, e
#   (B) as RECOMENDACOES que NAO foram aplicadas (dependem de verificacao sua),
#   (C) as CONFIRMACOES (coisas que pareciam problema mas estao corretas).
#
# Criterio adotado: so alterei o que tenho CERTEZA que e seguro (nao muda
# nenhum resultado ja publicado). Tudo que poderia mudar numero/figura ficou
# como recomendacao para VOCE decidir.
# =============================================================================


# =============================================================================
# 0. CONTEXTO IMPORTANTE (define as prioridades)
# =============================================================================
# O arquivo reference/expected_hashes.csv (a "prova" de reproducao byte-a-byte)
# cobre APENAS:
#     db/normalized_dim_location.csv
#     db/normalized_fact_cases.csv
#     db/normalized_fact_climate.csv
#     db/normalized_fact_population.csv
#     clinical_tables/clinical_demographic_long.csv
#     clinical_tables/clinical_{DENV,CHIKV,ZIKV}.xlsx
#
# >> Ou seja: a byte-identidade e exigida no ETL (etl/ + build_source_data) e
#    nas tabelas clinicas. Os CSVs de outputs/ (Mann-Kendall, Spearman, GLM,
#    Moran, bootstrap, etc.) NAO sao verificados por hash. Por isso, divergencias
#    numericas nos scripts de analise importam para a REPRODUTIBILIDADE
#    CIENTIFICA, mas nao quebram o teste de hash.


# =============================================================================
# 1. ALTERACOES APLICADAS  (seguras - ja gravadas na sua pasta)
# =============================================================================

# -----------------------------------------------------------------------------
# R/utils.R  ->  funcao kahan_by()   - APROVADO POR SOFIA
# -----------------------------------------------------------------------------
# PROBLEMA: o argumento default G = max(g) e um "promise" que o R so avalia no
# primeiro uso de G - e isso acontecia DEPOIS da linha `g <- g[ok]` (filtro de
# NA). Entao max(g) enxergava o g ja filtrado. Se o(s) grupo(s) de maior indice
# tivessem TODOS os valores NA, eles sumiriam do resultado (vetor curto demais).
#
# POR QUE E SEGURO: os 3 locais que chamam kahan_by (etl/build_climate.R, linhas
# ~112, ~113, ~128) JA passam o G explicitamente. Logo, a correcao e um NO-OP
# para o pipeline atual - NAO muda nenhum CSV nem os hashes. E so blindagem.
#
# ANTES:
#   kahan_by <- function(x, g, G = max(g), fun = c("sum", "mean")) {
#     fun <- match.arg(fun)
#     ok <- !is.na(x)
#
# DEPOIS:
#   kahan_by <- function(x, g, G = max(g), fun = c("sum", "mean")) {
#     fun <- match.arg(fun)
#     force(G)   # avalia max(g) com o g ORIGINAL, antes do filtro de NA abaixo
#     ok <- !is.na(x)


# =============================================================================
# 2. RECOMENDACOES  (NAO aplicadas - decisao sua; podem mudar resultado/figura)
# =============================================================================
# Severidade: [RISCO] = pode divergir do Python ou travar em certos dados
#             [MENOR] = borda rara / cosmetico / documentacao

# ---- [RISCO] -----------------------------------------------------------------

# R1. pipelines/build_source_data.R  (~linha 236)   conf: media
#     "Classfication" (CHIKV) escrito diferente de "Classification" (ZIKV).
#     find_table() usa startsWith() e da stop() se nao achar -> pode TRAVAR se o
#     titulo real da tabela no .pzfx estiver escrito certo.
#       list("Chikungunya", "A", c("Classfication", "Criteria", "Evolution")),
#     ACAO: conferir o titulo real da tabela no arquivo .pzfx do Chikungunya;
#           se for "Classification", corrigir a grafia. --------------------------- ACHO BOM VERIFICAR 

# R2. R/stats.R  (~linha 62)   natural_breaks()   conf: media
#     Usa style = "jenks" do classInt. O comentario cita mapclassify
#     NaturalBreaks E FisherJenks - mas no Python sao algoritmos DIFERENTES:
#     FisherJenks e otimo-exato (equivalente no R = style = "fisher"); o
#     NaturalBreaks e iterativo (nao garante otimo). "jenks" e "fisher" podem
#     dar quebras diferentes -> muda as classes de cor da Figura 6 (mapa
#     combinado) e das quebras geoespaciais.
#       b <- classInt::classIntervals(vals, n = k, style = "jenks")$brks
#     ACAO: conferir qual funcao o Python usou. Se FisherJenks, trocar para
#           style = "fisher".-------------------------------------------------- OLHAR EM PYTHON

# R3. pipelines/glm_common.R  (~linha 172)  fit_negbin()   conf: media
#     Estima o alpha do Negative Binomial com stats::optimize (Brent). O scipy
#     usa minimize_scalar(method="bounded"), outra implementacao -> o alpha pode
#     diferir nos ultimos digitos, e isso se propaga para coef/SE/IRR/CI do
#     NegBin (glm_negbin_*.csv). Nao e hash-checked, mas afeta reproducao exata.
#       opt <- suppressWarnings(stats::optimize(neg_llf, c(lo, hi), tol = 1e-6))
#     ACAO: conferir se o alpha bate (uns 8-10 digitos) com o resultado Python.------------- DESNECESSARIO

# R4. pipelines/glm_common.R  (~linha 226)  summarize_fit()   conf: baixa
#     Significant = unname(p) < ALPHA  ->  vira NA (celula vazia) quando p e NA
#     (ajuste degenerado / SE zero). No pandas, p < 0.05 daria False.
#     SUGESTAO: Significant = !is.na(p) & p < ALPHA
#     (so muda algo se existir p NA nos seus dados - verificar se e alcancavel.)------------- DESNECESSARIO

# R5. pipelines/bootstrap_ci.R  (~linha 44)  ci_row()   conf: baixa
#     se <- stats::sd(x) / sqrt(length(x))   ->  sd() usa ddof = 1.
#     As colunas Normal_CI_* sao DETERMINISTICAS (deveriam bater exatamente). Se
#     o Python calculou o SE com np.std(ddof = 0), diverge em toda linha.
#     ACAO: conferir o ddof no Python; se numpy default, usar
#           sqrt(mean((x - mean(x))^2)) / sqrt(length(x)).------------------ OLHAR EM PYTHON

# R6. pipelines/graphpad_pipeline.R  (~linha 386)   conf: media
#     writeLines(..., "figure_exclusions_captions.txt", useBytes = TRUE) grava
#     com LF (sep padrao). write_csv_py grava CRLF "as run on Windows"; se o
#     Python gravou este .txt em modo texto no Windows, era CRLF -> difere.
#     (Nao esta no expected_hashes.csv -> prioridade baixa.)
#     SUGESTAO: writeLines(..., sep = "\r\n") se quiser bater byte-a-byte.------------- DESNECESSARIO

# R7. pipelines/glm_common.R  (~linha 96)  wls_solve()   conf: baixa
#     Resolve por QR (qr.coef); statsmodels WLS usa pinv (SVD). Em designs quase
#     colineares podem diferir no ultimo ULP (exposto pelo py_float_str).
#     ACAO: so agir se observar divergencia; ai usar svd/pinv.----------------------- DESNECESSARIO

# R8. pipelines/additional/adv_dlnm.R  (~linha 98)   conf: media
#     V <- solve(crossprod(X, X * oim_factor(...)))  sem tryCatch -> TRAVA
#     ("system is computationally singular") se X'WX for singular; statsmodels
#     usa pseudo-inversa. SUGESTAO:
#       XtWX <- crossprod(X, X * oim_factor(sub$Cases, f$mu, nb$alpha))
#       V <- tryCatch(solve(XtWX), error = function(e) MASS::ginv(XtWX))--------------------------- ACHO BOM VERIFICAR 

# R9. etl/build_cases_clinical.R  (~linha 64)  is_py_digit()   conf: baixa
#     grepl("^[0-9...]+$", v, perl = TRUE) com padrao UTF-8 sobre strings
#     marcadas como Latin-1 -> borda para codigos com superscript. So afeta
#     rotulo de categoria raro.  SUGESTAO: usar enc2utf8(v) antes do grepl.----------------------- DESNECESSARIO

# R10. etl/build_climate.R  (~linha 34)  parse_decimal_comma()   conf: baixa
#     Le INMET com strip.white = FALSE; tokens com espaco (" 21,5") viram NA,
#     enquanto o parser do pandas apara espaco -> " 21,5" = 21.5. Se os arquivos
#     tiverem campos com espaco, Temp/Rain divergem. SUGESTAO: trimws() antes.--------------------------- ACHO BOM VERIFICAR 

# ---- [MENOR] -----------------------------------------------------------------

# M1. bootstrap_ci.R: a Media usa x[!is.na(x)] (mantem Inf) e o bootstrap usa
#     is.finite (remove Inf) -> inconsistencia latente se houver Inf. Filtrar
#     ambos com is.finite.---------------------------------------------------------- DESNECESSARIO
# M2. additional/adv_cocirculation.R (~86): o df do teste-t usa o n inteiro de
#     "a", ignorando o shift do lag e os NA de "b" -> p enviesado p/ baixo.
#     (Documentado como "naive", mas o n esta errado p/ o que e testado.)----------------- ACHO BOM VERIFICAR 
# M3. additional/adv_forecast.R: o cabecalho diz "forward-difference / so
#     L-BFGS-B", mas o codigo usa diferenca CENTRAL + passo BFGS. So comentario.---------- DESNECESSARIO
# M4. additional/adv_covid_its.R (~38): t_since_pandemic so e != 0 DENTRO da
#     janela (pulso), nao e uma tendencia pos-intervencao permanente. Nome/doc.---------------- ACHO BOM VERIFICAR 
# M5. additional/adv_common.R (~13): RUN_DATE <- "2026-10-08" fixo -> sempre
#     grava em outputs/2026-10-08/. Confirmar se e proposital (vs Sys.Date()).--------- DESNECESSARIO
# M6. etl/build_climate.R (~141): lags calculados sobre valores JA arredondados
#     (o py_round externo vira no-op). Provavelmente correto (pandas arredonda
#     antes), mas confirmar a ordem no Python (bordas .xx5).--------- ACHO BOM VERIFICAR 
# M7. R/utils.R py_float_str(): aproxima o repr() de float do Python com a
#     escada 15/16/17g; diverge em magnitudes extremas / fronteira 1e16. OK para
#     taxas, p-valores e incidencias (validado pelo expected_hashes.csv).--------- ACHO BOM VERIFICAR 
# M8. R/utils.R py_strip(): cobre espaco/tab/NBSP/separadores/NEL, mas nao os
#     espacos Unicode exoticos (\u2000-\u200a, \u3000). Improvavel nos seus dados.--------- DESNECESSARIO
# M9. R/utils.R write_xlsx_sheets(): corta nome de aba em 31 chars; se dois
#     nomes ficarem iguais apos o corte, o Excel reclama. Conferir nomes longos.--------- DESNECESSARIO
# M10. pipelines/arbovirus_geospatial.R: usa caracteres nao-ASCII (- longo e ->)
#     no texto dos graficos -> salvar o arquivo em UTF-8 para nao quebrar no
#     Windows. Alem disso, "moran_global" e usado como funcao E como variavel
#     (sombra) em main()/plot_combined_map() - funciona, mas e confuso.--------- DESNECESSARIO

