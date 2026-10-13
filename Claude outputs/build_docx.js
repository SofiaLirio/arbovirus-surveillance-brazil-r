const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, BorderStyle, ShadingType,
} = require("docx");

const BRAND = "2E5A87";
const MONO = "Consolas";

function h1(t) { return new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 280, after: 120 }, children: [new TextRun({ text: t, bold: true, color: BRAND })] }); }
function h2(t) { return new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 200, after: 80 }, children: [new TextRun({ text: t, bold: true })] }); }
function p(runs) { return new Paragraph({ spacing: { after: 100 }, children: Array.isArray(runs) ? runs : [new TextRun(runs)] }); }
function t(text, opts = {}) { return new TextRun({ text, ...opts }); }

// code block: each line a monospace paragraph on light-grey shading
function code(lines) {
  const arr = Array.isArray(lines) ? lines : lines.split("\n");
  return arr.map((ln, i) => new Paragraph({
    spacing: { after: i === arr.length - 1 ? 100 : 0, before: i === 0 ? 40 : 0 },
    shading: { type: ShadingType.CLEAR, fill: "F2F2F2" },
    children: [new TextRun({ text: ln === "" ? " " : ln, font: MONO, size: 17 })],
  }));
}

function bullet(text) { return new Paragraph({ bullet: { level: 0 }, spacing: { after: 60 }, children: [new TextRun(text)] }); }

// ---- recommendations table ----
const NB = { style: BorderStyle.SINGLE, size: 2, color: "BBBBBB" };
const BORDERS = { top: NB, bottom: NB, left: NB, right: NB, insideHorizontal: NB, insideVertical: NB };
const COLW = [560, 1700, 820, 3380, 2300, 600]; // sums 9360
const TABLEW = COLW.reduce((a, b) => a + b, 0);

function cell(text, w, opts = {}) {
  return new TableCell({
    width: { size: w, type: WidthType.DXA },
    shading: opts.head ? { type: ShadingType.CLEAR, fill: BRAND } : undefined,
    margins: { top: 40, bottom: 40, left: 70, right: 70 },
    children: (Array.isArray(text) ? text : [text]).map(s =>
      new Paragraph({ children: [new TextRun({ text: s, bold: !!opts.head, color: opts.head ? "FFFFFF" : (opts.color || "000000"), size: opts.head ? 17 : 16 })] })),
  });
}
function row(cells) { return new TableRow({ children: cells }); }
function headRow() {
  const hs = ["ID", "Arquivo (linha)", "Sev.", "Problema", "Sugestão", "Conf."];
  return new TableRow({ tableHeader: true, children: hs.map((hdr, i) => cell(hdr, COLW[i], { head: true })) });
}
function recRow(r) {
  const sevColor = r.sev.startsWith("RISCO") ? "B00020" : "7A5C00";
  return row([
    cell(r.id, COLW[0]),
    cell(r.file, COLW[1]),
    cell(r.sev, COLW[2], { color: sevColor }),
    cell(r.prob, COLW[3]),
    cell(r.fix, COLW[4]),
    cell(r.conf, COLW[5]),
  ]);
}

const RECS = [
  { id: "R1", file: "build_source_data.R (~236)", sev: "RISCO", prob: "\"Classfication\" (CHIKV) difere de \"Classification\" (ZIKV); find_table() da stop() se nao achar -> pode travar.", fix: "Conferir titulo real no .pzfx; corrigir grafia se preciso.", conf: "media" },
  { id: "R2", file: "R/stats.R (~62) natural_breaks", sev: "RISCO", prob: "style=\"jenks\"; se o Python usou FisherJenks, o equivalente exato e \"fisher\". Muda classes de cor (Fig. 6/mapas).", fix: "Conferir funcao do Python; se FisherJenks, usar style=\"fisher\".", conf: "media" },
  { id: "R3", file: "glm_common.R (~172) fit_negbin", sev: "RISCO", prob: "alpha do NegBin via stats::optimize (Brent) vs scipy minimize_scalar bounded -> alpha pode diferir nos ultimos digitos (propaga p/ coef/IRR/CI).", fix: "Conferir se o alpha bate ~8-10 digitos com o Python.", conf: "media" },
  { id: "R4", file: "glm_common.R (~226)", sev: "RISCO", prob: "Significant = p < ALPHA vira NA quando p e NA; pandas daria False.", fix: "Significant = !is.na(p) & p < ALPHA (se alcancavel).", conf: "baixa" },
  { id: "R5", file: "bootstrap_ci.R (~44)", sev: "RISCO", prob: "se = sd(x)/sqrt(n) usa ddof=1; Normal_CI e deterministico. Se o Python usou np.std(ddof=0), diverge.", fix: "Conferir ddof; se numpy, usar sqrt(mean((x-mean)^2))/sqrt(n).", conf: "baixa" },
  { id: "R6", file: "graphpad_pipeline.R (~386)", sev: "RISCO", prob: "figure_exclusions_captions.txt gravado com LF; Windows texto seria CRLF. (Nao esta no expected_hashes.)", fix: "writeLines(..., sep=\"\\r\\n\") se quiser byte-identico.", conf: "media" },
  { id: "R7", file: "glm_common.R (~96) wls_solve", sev: "RISCO", prob: "Resolve por QR; statsmodels WLS usa pinv (SVD). Diferenca de ULP em designs quase colineares.", fix: "So agir se houver divergencia; usar svd/pinv.", conf: "baixa" },
  { id: "R8", file: "additional/adv_dlnm.R (~98)", sev: "RISCO", prob: "V <- solve(...) sem tryCatch -> trava se X'WX singular; statsmodels usa pseudo-inversa.", fix: "V <- tryCatch(solve(XtWX), error=\\(e) MASS::ginv(XtWX)).", conf: "media" },
  { id: "R9", file: "build_cases_clinical.R (~64)", sev: "RISCO", prob: "is_py_digit: padrao UTF-8 sobre strings Latin-1 -> borda p/ codigos com superscript (rotulo raro).", fix: "grepl(..., enc2utf8(v), perl=TRUE).", conf: "baixa" },
  { id: "R10", file: "build_climate.R (~34)", sev: "RISCO", prob: "parse_decimal_comma com strip.white=FALSE: \" 21,5\" vira NA; pandas apara espaco.", fix: "trimws() antes do regex.", conf: "baixa" },
  { id: "M1", file: "bootstrap_ci.R (~33/43)", sev: "MENOR", prob: "Media mantem Inf (!is.na), bootstrap remove Inf (is.finite) -> inconsistencia se houver Inf.", fix: "Filtrar ambos com is.finite.", conf: "baixa" },
  { id: "M2", file: "additional/adv_cocirculation.R (~86)", sev: "MENOR", prob: "df do teste-t usa n inteiro de 'a', ignora shift do lag e NA de 'b' -> p enviesado p/ baixo.", fix: "Usar o n do overlap real no best_lag.", conf: "media" },
  { id: "M3", file: "additional/adv_forecast.R", sev: "MENOR", prob: "Docstring diz forward-difference/L-BFGS-B; codigo usa diferenca central + BFGS. So comentario.", fix: "Alinhar o comentario ao codigo.", conf: "alta" },
  { id: "M4", file: "additional/adv_covid_its.R (~38)", sev: "MENOR", prob: "t_since_pandemic so !=0 dentro da janela (pulso), nao tendencia permanente. Risco de leitura.", fix: "Renomear coluna ou documentar interpretacao.", conf: "baixa" },
  { id: "M5", file: "additional/adv_common.R (~13)", sev: "MENOR", prob: "RUN_DATE fixo \"2026-10-08\" -> sempre grava na mesma pasta.", fix: "Confirmar se e proposital (vs Sys.Date()).", conf: "alta" },
  { id: "M6", file: "build_climate.R (~141)", sev: "MENOR", prob: "Lags calculados sobre valores ja arredondados. Provavelmente correto.", fix: "Confirmar ordem (arredonda-entao-shift) no Python.", conf: "media" },
  { id: "M7", file: "R/utils.R py_float_str", sev: "MENOR", prob: "Aproxima repr() do Python (15/16/17g); diverge em magnitudes extremas / fronteira 1e16.", fix: "OK p/ taxas e p-valores (validado por hashes).", conf: "media" },
  { id: "M8", file: "R/utils.R py_strip", sev: "MENOR", prob: "Nao cobre espacos Unicode exoticos (\\u2000-\\u200a, \\u3000).", fix: "Improvavel nos dados; so se aparecer.", conf: "baixa" },
  { id: "M9", file: "R/utils.R write_xlsx_sheets", sev: "MENOR", prob: "Corte de nome de aba em 31 chars pode colidir -> Excel reclama.", fix: "Conferir nomes de aba muito longos/parecidos.", conf: "baixa" },
  { id: "M10", file: "arbovirus_geospatial.R", sev: "MENOR", prob: "Caracteres nao-ASCII nos graficos (salvar UTF-8); moran_global usado como funcao e variavel (sombra).", fix: "Salvar em UTF-8; renomear a variavel.", conf: "media" },
];

const doc = new Document({
  creator: "Claude (Cowork)",
  title: "Relatorio de Revisao - arbovirus-surveillance-brazil-r",
  styles: { default: { document: { run: { font: "Calibri", size: 21 } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    children: [
      new Paragraph({ spacing: { after: 60 }, children: [new TextRun({ text: "Relatório de Revisão de Scripts", bold: true, size: 40, color: BRAND })] }),
      new Paragraph({ spacing: { after: 40 }, children: [new TextRun({ text: "Projeto: arbovirus-surveillance-brazil-r", size: 24, bold: true })] }),
      p([t("Revisão por: Claude (Cowork)   •   Solicitado por: Sofia Lirio (LII LAB)   •   Data: 2026-10-13")]),
      new Paragraph({ border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: BRAND } }, spacing: { after: 160 }, children: [] }),

      h1("Resumo executivo"),
      p([t("Foram verificados "), t("36 scripts R", { bold: true }), t(" do projeto (R/, pipelines/, pipelines/additional/, etl/, tools/). O código é, no geral, um port muito cuidadoso de um pipeline Python (pandas/numpy/statsmodels) para epidemiologia de arboviroses, com meta de reprodução byte-a-byte dos CSVs.")]),
      bullet("2 alterações aplicadas (seguras — não mudam resultados já publicados)."),
      bullet("20 recomendações para você verificar (NÃO aplicadas) — ver tabela."),
      bullet("6 confirmações (pontos que pareciam problema mas estão corretos)."),
      p([t("Critério: só alterei o que tenho certeza que é seguro. Tudo que poderia mudar um número ou figura ficou como recomendação, para a sua decisão.", { italics: true })]),

      h1("0. Contexto que define as prioridades"),
      p([t("O "), t("reference/expected_hashes.csv", { font: MONO, size: 18 }), t(" (a prova de reprodução byte-a-byte) cobre APENAS os CSVs de "), t("db/", { font: MONO, size: 18 }), t(" e as tabelas clínicas. Os CSVs de "), t("outputs/", { font: MONO, size: 18 }), t(" (Mann-Kendall, Spearman, GLM, Moran, bootstrap) NÃO são verificados por hash.")]),
      p([t("Consequência: a byte-identidade é crítica no ETL ("), t("etl/", { font: MONO, size: 18 }), t(" + "), t("build_source_data.R", { font: MONO, size: 18 }), t("). Divergências nos scripts de análise afetam a reprodutibilidade científica, mas não quebram o teste de hash.")]),

      h1("1. Alterações aplicadas (já gravadas na sua pasta)"),

      h2("[A1] R/utils.R — kahan_by(): reforço defensivo (NO-OP)"),
      p([t("Problema: o default "), t("G = max(g)", { font: MONO, size: 18 }), t(" só é avaliado no 1º uso de G, que acontecia DEPOIS do filtro de NA ("), t("g <- g[ok]", { font: MONO, size: 18 }), t("). Se o grupo de maior índice tivesse todos os valores NA, ele sumiria do resultado.")]),
      p([t("Por que é seguro: os 3 locais que chamam kahan_by (em "), t("etl/build_climate.R", { font: MONO, size: 18 }), t(") já passam o G explicitamente. Logo a correção é um "), t("no-op", { bold: true }), t(" — não muda nenhum CSV nem os hashes; é só blindagem.")]),
      p([t("Antes:", { bold: true })]),
      ...code(["kahan_by <- function(x, g, G = max(g), fun = c(\"sum\", \"mean\")) {", "  fun <- match.arg(fun)", "  ok <- !is.na(x)"]),
      p([t("Depois:", { bold: true })]),
      ...code(["kahan_by <- function(x, g, G = max(g), fun = c(\"sum\", \"mean\")) {", "  fun <- match.arg(fun)", "  force(G)   # avalia max(g) com o g ORIGINAL, antes do filtro de NA", "  ok <- !is.na(x)"]),

      h2("[A2] R/AnotaçõesSofia.R — remoção de \"-\" solto (linha 170)"),
      p([t("Problema: um "), t("-", { font: MONO, size: 18 }), t(" sozinho numa linha (resquício de colagem) entre "), t("order_c()", { font: MONO, size: 18 }), t(" e "), t("py_float_str()", { font: MONO, size: 18 }), t(" quebrava o source() do arquivo. (É seu arquivo de estudo, não usado pelo pipeline.)")]),
      p([t("Antes:", { bold: true })]),
      ...code(["order_c <- function(...) order(..., method = \"radix\")", "-", "", "py_float_str <- function(x) {"]),
      p([t("Depois:", { bold: true })]),
      ...code(["order_c <- function(...) order(..., method = \"radix\")", "", "py_float_str <- function(x) {"]),

      h1("2. Recomendações (NÃO aplicadas — decisão sua)"),
      p([t("Severidade: "), t("RISCO", { bold: true, color: "B00020" }), t(" = pode divergir do Python ou travar em certos dados; "), t("MENOR", { bold: true, color: "7A5C00" }), t(" = borda rara / cosmético / documentação.")]),
      new Table({ width: { size: TABLEW, type: WidthType.DXA }, columnWidths: COLW, borders: BORDERS, rows: [headRow(), ...RECS.map(recRow)] }),

      h1("3. Confirmações (verificado — sem ação)"),
      bullet("Dependências do geospatial OK: data_loader.R (linha 19) faz source(R/stats.R) e (linha 142) define knn_weights(). Por isso moran_global/moran_local/knn_weights são encontrados. Não era bug."),
      bullet("R/stats.R: Mann-Kendall (correção de continuidade), Sen, Spearman (t do scipy), folded_p (= esda), Moran global/local (desvio-padrão populacional, permutação condicional) — matemática fiel."),
      bullet("R/utils.R: kahan_by reproduz a soma de Kahan do pandas; write_csv_py/py_col_str/csv_quote reproduzem o to_csv do Windows (CRLF, UTF-8, True/False, int sem .0); py_round = numpy round (bankers)."),
      bullet("arbovirus_pipeline_en.R: no Mann-Kendall, Significant = res$p < ALPHA é seguro (p sempre finito); o Spearman já usa !is.na(). O problema de NA (R4) não ocorre aqui."),
      bullet("ETL: age_band = pd.cut(right=FALSE), merges sem linhas órfãs, order_c = ordenação do pandas, verificação sha256 em check_db.R/compare_with_archive.R corretas."),
      bullet("additional/: epidemic_threshold, spatial_regression, forecast (Kalman/KFAS) e dlnm (crossbasis) conferidos e coerentes (fora R8)."),

      h1("4. Como carregar / validar"),
      p([t("Carregar só as funções (utils + stats) para explorar no Positron:")]),
      ...code([".root <- \"C:/Users/User/OneDrive/EBMSP/GRUPO VIVIANE/DENGUE/arbovirus-surveillance-brazil-r\"", "source(file.path(.root, \"R\", \"utils.R\"))", "source(file.path(.root, \"R\", \"stats.R\"))"]),
      p([t("Validar a byte-identidade do DB (o que o projeto garante):")]),
      ...code(["Rscript etl/check_db.R", "Rscript tools/compare_with_archive.R"]),
      p([t("Rodar um pipeline de análise inteiro (precisa do db/ e dos dados):")]),
      ...code(["Rscript pipelines/arbovirus_pipeline_en.R", "Rscript pipelines/arbovirus_geospatial.R"]),
    ],
  }],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync("/mnt/user-data/outputs/RELATORIO_REVISAO_SOFIA.docx", buf);
  console.log("ok, bytes:", buf.length);
});
