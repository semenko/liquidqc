#!/usr/bin/env Rscript
# Create synthetic test data for dupRust integration tests
# This script:
# 1. Creates a small GTF annotation file
# 2. Creates a synthetic BAM file with known duplicate patterns
# 3. Runs the original dupRadar to produce reference outputs

library(dupRadar)

cat("=== Creating test data for dupRust ===\n")

outdir <- "tests/data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create("tests/expected", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Create a synthetic GTF file
# ============================================================
cat("Creating GTF file...\n")

# Define genes on chr1 with various expression levels
# gene_id, chr, start, end, strand, exon_starts, exon_ends
genes <- list(
  # High expression gene - 3 exons, ~2000bp total
  list(id = "GENE_HIGH", chr = "chr1", strand = "+",
       exons = data.frame(start = c(1000, 2000, 3000),
                          end   = c(1500, 2500, 3800))),
  # Medium expression gene - 2 exons, ~1000bp total
  list(id = "GENE_MED", chr = "chr1", strand = "+",
       exons = data.frame(start = c(5000, 6000),
                          end   = c(5600, 6500))),
  # Low expression gene - 1 exon, ~500bp
  list(id = "GENE_LOW", chr = "chr1", strand = "-",
       exons = data.frame(start = c(8000),
                          end   = c(8500))),
  # Very low expression gene - 1 exon, ~300bp
  list(id = "GENE_VLOW", chr = "chr1", strand = "-",
       exons = data.frame(start = c(10000),
                          end   = c(10300))),
  # Zero expression gene (no reads)
  list(id = "GENE_ZERO", chr = "chr1", strand = "+",
       exons = data.frame(start = c(12000),
                          end   = c(12400))),
  # Another medium gene on chr2
  list(id = "GENE_MED2", chr = "chr2", strand = "+",
       exons = data.frame(start = c(1000, 2000),
                          end   = c(1400, 2600))),
  # Gene with overlapping exons (tests length calculation)
  list(id = "GENE_OVERLAP", chr = "chr2", strand = "-",
       exons = data.frame(start = c(4000, 4200, 4500),
                          end   = c(4300, 4600, 4800))),
  # Small gene for high RPK
  list(id = "GENE_SMALL", chr = "chr2", strand = "+",
       exons = data.frame(start = c(6000),
                          end   = c(6100)))
)

# Write GTF file
gtf_file <- file.path(outdir, "test.gtf")
gtf_lines <- c()
for (gene in genes) {
  # Gene line
  gene_start <- min(gene$exons$start)
  gene_end <- max(gene$exons$end)
  gtf_lines <- c(gtf_lines, sprintf(
    '%s\ttest\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_name "%s";',
    gene$chr, gene_start, gene_end, gene$strand, gene$id, gene$id))

  # Exon lines
  for (i in 1:nrow(gene$exons)) {
    gtf_lines <- c(gtf_lines, sprintf(
      '%s\ttest\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_name "%s"; exon_number "%d";',
      gene$chr, gene$exons$start[i], gene$exons$end[i], gene$strand,
      gene$id, gene$id, i))
  }
}
writeLines(gtf_lines, gtf_file)
cat("  GTF written to:", gtf_file, "\n")

# ============================================================
# 2. Create synthetic BAM file using SAM format
# ============================================================
cat("Creating BAM file...\n")

sam_file <- file.path(outdir, "test.sam")
bam_file <- file.path(outdir, "test.bam")

# SAM header
header <- c(
  "@HD\tVN:1.6\tSO:coordinate",
  "@SQ\tSN:chr1\tLN:20000",
  "@SQ\tSN:chr2\tLN:20000",
  "@RG\tID:sample1\tSM:sample1"
)

set.seed(42)

# Generate reads for each gene
# Columns: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, tags
reads <- c()
read_id <- 1

make_read <- function(name, flag, chr, pos, mapq = 255, cigar = "50M",
                      nh = 1, is_dup = FALSE) {
  if (is_dup) flag <- bitwOr(flag, 1024)
  seq <- paste(rep("A", 50), collapse = "")
  qual <- paste(rep("I", 50), collapse = "")
  tags <- sprintf("NH:i:%d\tHI:i:1\tRG:Z:sample1", nh)
  sprintf("%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\t%s",
          name, flag, chr, pos, mapq, cigar, seq, qual, tags)
}

# GENE_HIGH: 200 total reads, 40 duplicates (20% dup rate)
for (i in 1:200) {
  is_dup <- i > 160
  pos <- sample(c(1000:1450, 2000:2450, 3000:3750), 1)
  reads <- c(reads, make_read(
    sprintf("read_%05d", read_id), 0, "chr1", pos, nh = 1, is_dup = is_dup))
  read_id <- read_id + 1
}

# GENE_MED: 50 total reads, 10 duplicates (20% dup rate)
for (i in 1:50) {
  is_dup <- i > 40
  pos <- sample(c(5000:5550, 6000:6450), 1)
  reads <- c(reads, make_read(
    sprintf("read_%05d", read_id), 0, "chr1", pos, nh = 1, is_dup = is_dup))
  read_id <- read_id + 1
}

# GENE_LOW: 10 total reads, 5 duplicates (50% dup rate - suspicious)
for (i in 1:10) {
  is_dup <- i > 5
  pos <- sample(8000:8450, 1)
  reads <- c(reads, make_read(
    sprintf("read_%05d", read_id), 0, "chr1", pos, nh = 1, is_dup = is_dup))
  read_id <- read_id + 1
}

# GENE_VLOW: 3 total reads, 2 duplicates (67% dup rate - very suspicious)
for (i in 1:3) {
  is_dup <- i > 1
  pos <- sample(10000:10250, 1)
  reads <- c(reads, make_read(
    sprintf("read_%05d", read_id), 0, "chr1", pos, nh = 1, is_dup = is_dup))
  read_id <- read_id + 1
}

# GENE_MED2: 80 total reads, 15 duplicates, plus 10 multimappers
for (i in 1:80) {
  is_dup <- i > 65
  pos <- sample(c(1000:1350, 2000:2550), 1)
  reads <- c(reads, make_read(
    sprintf("read_%05d", read_id), 0, "chr2", pos, nh = 1, is_dup = is_dup))
  read_id <- read_id + 1
}
# Multimappers for GENE_MED2
for (i in 1:10) {
  pos <- sample(c(1000:1350, 2000:2550), 1)
  reads <- c(reads, make_read(
    sprintf("read_%05d", read_id), 0, "chr2", pos, nh = 3, is_dup = FALSE))
  read_id <- read_id + 1
}

# GENE_OVERLAP: 30 total reads, 5 duplicates
for (i in 1:30) {
  is_dup <- i > 25
  pos <- sample(4000:4750, 1)
  reads <- c(reads, make_read(
    sprintf("read_%05d", read_id), 0, "chr2", pos, nh = 1, is_dup = is_dup))
  read_id <- read_id + 1
}

# GENE_SMALL: 100 total reads, 60 duplicates (60% dup rate - expected, small gene)
for (i in 1:100) {
  is_dup <- i > 40
  pos <- sample(6000:6050, 1)
  reads <- c(reads, make_read(
    sprintf("read_%05d", read_id), 0, "chr2", pos, nh = 1, is_dup = is_dup))
  read_id <- read_id + 1
}

# Also add some unmapped reads (should be excluded from N)
for (i in 1:5) {
  reads <- c(reads, make_read(
    sprintf("unmapped_%05d", read_id), 4, "*", 0, mapq = 0, nh = 0))
  read_id <- read_id + 1
}

# Write SAM file
writeLines(c(header, reads), sam_file)

# Convert SAM to sorted BAM
system2("samtools", c("view", "-bS", sam_file, "-o",
                       file.path(outdir, "test_unsorted.bam")))
system2("samtools", c("sort", file.path(outdir, "test_unsorted.bam"),
                       "-o", bam_file))
system2("samtools", c("index", bam_file))
file.remove(sam_file, file.path(outdir, "test_unsorted.bam"))
cat("  BAM written to:", bam_file, "\n")

# ============================================================
# 3. Run dupRadar to produce reference outputs
# ============================================================
cat("Running dupRadar...\n")

dm <- analyzeDuprates(
  bam = bam_file,
  gtf = gtf_file,
  stranded = 0,
  paired = FALSE,
  threads = 1,
  verbose = TRUE
)

cat("\nDuplication matrix:\n")
print(dm)
cat("\n")

# Save duplication matrix
dm_file <- "tests/expected/dupMatrix.txt"
write.table(dm, file = dm_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Matrix written to:", dm_file, "\n")

# Fit model
fit <- duprateExpFit(dm)
cat("\nFit results:\n")
cat("  intercept:", fit$intercept, "\n")
cat("  slope:", fit$slope, "\n")
cat("  beta0:", coef(fit$glm)[1], "\n")
cat("  beta1:", coef(fit$glm)[2], "\n")

# Save fit
fit_file <- "tests/expected/intercept_slope.txt"
writeLines(c(
  paste("intercept", fit$intercept, sep = "\t"),
  paste("slope", fit$slope, sep = "\t")
), fit_file)
cat("  Fit written to:", fit_file, "\n")

# Save stats
stats <- getDupMatStats(dm)
cat("\nStats:\n")
print(stats)
stats_file <- "tests/expected/stats.txt"
write.table(t(as.data.frame(stats)), file = stats_file, sep = "\t",
            row.names = TRUE, col.names = FALSE, quote = FALSE)

# Generate plots
pdf("tests/expected/duprateExpDens.pdf", width = 8, height = 6)
duprateExpDensPlot(dm)
dev.off()
cat("  Density plot written\n")

pdf("tests/expected/duprateExpBoxplot.pdf", width = 8, height = 6)
duprateExpBoxplot(dm)
dev.off()
cat("  Boxplot written\n")

pdf("tests/expected/expressionHist.pdf", width = 8, height = 6)
expressionHist(dm)
dev.off()
cat("  Expression histogram written\n")

# Save fitted curve data for comparison
x_vals <- seq(min(log10(dm$RPK[dm$RPK > 0])),
              max(log10(dm$RPK[dm$RPK > 0])), length.out = 1000)
y_vals <- predict(fit$glm, data.frame(x = x_vals), type = "response")
curve_data <- data.frame(log10RPK = x_vals, dupRate = y_vals)
write.table(curve_data, file = "tests/expected/fit_curve.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Fit curve data written\n")

cat("\n=== Test data generation complete ===\n")
