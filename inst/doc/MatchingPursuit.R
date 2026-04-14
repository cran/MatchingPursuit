## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE)

## ----old_par, include = FALSE-------------------------------------------------
old.par <- par("mfrow", "mai", "pty")

## ----R_version, include = TRUE, echo = FALSE----------------------------------
out <- version
cat(
  "R version:    ", 
  out$major, ".", 
  out$minor, 
  "\n", 
  "Generated on: ", 
  format(Sys.time(), "%d-%B-%Y"), 
  sep = ""
)

## ----load_library-------------------------------------------------------------
library(MatchingPursuit)

## ----empi_install, include = TRUE, echo = TRUE--------------------------------
# empi.install()

## ----empi_find, include = TRUE, echo = TRUE-----------------------------------
empi.check()

## ----7_non_stationary, include = TRUE, echo = TRUE, warning = FALSE-----------
fs <- 1024
T <- 1
t <- seq(0, T - 1 / fs, 1 / fs)
N <- length(t)

# 7 non-stationary signals.
x1 <- sin(2 * pi * (10 + 40 * t) * t)                            # linear chirp
x2 <- sin(2 * pi * (20 * t^2) * t)                               # nonlinear chirp
x3 <- (1 + 0.5 * sin(2 * pi * 2 * t))  *  sin(2 * pi * 30 * t)   # AM
x4 <- sin(2 * pi * 50 * t + 5 * sin(2 * pi * 3 * t))             # FM
x5 <- exp(-2 * t)  *  sin(2 * pi * 60 * t)                       # decreasing amplitude
x6 <- sin(2 * pi * (5 + 20 * sin(2 * pi * t)) * t)               # frequency modulated sine wave
x7 <- t * sin(2 * pi * 40 * t)                                   # increasing amplitude

signal <- data.frame(x = x1 + x2 + x3 + x4 + x5 + x6 + x7)

## ----7_non_stationary_plot, include = TRUE, echo = FALSE, warning = FALSE, fig.width = 7, fig.height = 7----
range <- range(signal)

par(mfrow = c(8, 1), pty = "m", mai = c(0.2, 0.4, 0.2, 0.1))

plot(t, signal$x, type = "l", col = "blue", main = "", xlab = "", ylab = "")
plot(t, x1, type = "l", ylab = "", xlab = "", xaxt = "n", yaxt = "s", main = "Linear chirp")
plot(t, x2, type = "l", ylab = "", xlab = "", xaxt = "n", yaxt = "s", main = "Nonlinear chirp")
plot(t, x3, type = "l", ylab = "", xlab = "", xaxt = "n", yaxt = "s", main = "AM")
plot(t, x4, type = "l", ylab = "", xlab = "", xaxt = "n", yaxt = "s", main = "FM")
plot(t, x5, type = "l", ylab = "", xlab = "", xaxt = "n", yaxt = "s", main = "Decreasing amplitude")
plot(t, x6, type = "l", ylab = "", xlab = "", xaxt = "n", yaxt = "s", main = "Frequency modulated sine wave")
plot(t, x7, type = "l", ylab = "", xlab = "", xaxt = "n", yaxt = "s", main = "Increasing amplitude")

par(old.par)

## ----signal, include = TRUE, echo = TRUE, warning = FALSE---------------------
file <- system.file("extdata", "sample1.csv", package = "MatchingPursuit")

out <- read.csv(file, header = FALSE)
head(out)

out <- read.csv.signals(file)

str(out)

head(out$signal)

out$sampling.rate

## ----empi_execute, include = TRUE, echo = TRUE, warning = FALSE---------------
# file <- system.file("extdata", "sample1.csv", package = "MatchingPursuit")
# signal <- read.csv.signals(file)

# empi.out <- empi.execute (
#   signal = signal,
#   empi.options = "-o local --gabor -i 25",
#   write.to.file = TRUE,
#   path = NULL,
#   file.name = "sample1.db"
# )

## ----sample1, include = TRUE, echo = TRUE, warning = FALSE, fig.width = 7, fig.height = 7----
# Reading a SQLite file in which all generated atom parameters are stored.
file <- system.file("extdata", "sample1.db", package = "MatchingPursuit")

# Create time-frequency map based on MP atoms.
out <- empi2tf(
  db.file = file,
  channel = 1,
  mode = "sqrt",
  freq.divide = 4,
  increase.factor = 4,
  display.crosses = FALSE,
  display.atom.numbers = TRUE,
  out.mode = "plot"
)

## ----eeg_edf_reading, include = TRUE, echo = TRUE, warning = FALSE------------
file <- system.file("extdata", "EEG.edf", package = "MatchingPursuit")

# Read signal parameters and display them in a tabular form.
read.edf.params(file)

## ----eeg_edf_reading_resampling, include = TRUE, echo = TRUE, warning = FALSE----
# Original signal
out1 <- read.edf.signals(file, resampling = FALSE, from = 0, to = 10)
signal <- out1$signals
sampling.rate <- out1$sampling.rate
sampling.rate

# Resampled signal
out2 <- read.edf.signals(file, resampling = TRUE, f.new = 64,  from = 0, to = 10)
signal.resampled <- out2$signals
sampling.rate.resampled <- out2$sampling.rate
sampling.rate.resampled

## ----sampling, include = TRUE, echo = FALSE, warning = FALSE, fig.width = 7, fig.height = 4----
par(mfrow = c(2, 1), pty = "m", mai = c(0.8, 0.5, 0.5, 0.5))

# Not-filtered signal (raw signal).
plot(
  signal[, 1],
  type = "l",
  panel.first = grid(),
  main = "Original signal (256 Hz), channel #1",
  xlab = "sample points",
  ylab = "",
  col = "blue"
)

# Signal after filtering.
plot(
  signal.resampled[, 1],
  type = "l",
  panel.first = grid(),
  main = "Signal after downsampling (64 Hz), channel #1",
  xlab = "sample points",
  ylab = "",
  col = "blue"
)

par(old.par)

## ----montage, include = TRUE, echo = TRUE, warning = FALSE--------------------
# Pairs of signals for bipolar montage (so called "double banana").
pairs <- list(
  c("Fp2", "F4"), c("F4", "C4"), c("C4", "P4"), c("P4", "O2"), c("Fp1", "F3"), c("F3", "C3"),  
  c("C3", "P3"), c("P3", "O1"), c("Fp2", "F8"), c("F8", "T4"), c("T4", "T6"), c("T6", "O2"),
  c("Fp1", "F7"), c("F7", "T3"), c("T3", "T5"), c("T5", "O1"), c("Fz", "Cz"), c("Cz", "Pz")
)

# Make the bipolar montage.
bip.montage <- eeg.montage(signal, montage.type = c("bipolar"), bipolar.pairs = pairs)

# Original signal.
head(signal[, 1:5])

# Signal after banana montage.
head(bip.montage[, 1:5])


## ----filtering, include = TRUE, echo = TRUE, warning = FALSE------------------
# Filter parameters that will be used (quite typical in filtering EEG signals).
fc <- filters.coeff(
   fs = sampling.rate,
   notch = c(49, 51),
   lowpass = 40,
   highpass = 1,
)

# Filtering input signals.
bip.montage.filt <- bip.montage 

for (m in 1:ncol(bip.montage)) {
  bip.montage.filt[, m] = signal::filtfilt(fc$notch, bip.montage[, m])         # 50Hz notch filter
  bip.montage.filt[, m] = signal::filtfilt(fc$lowpass, bip.montage.filt[, m])  # Low pass IIR Butterworth
  bip.montage.filt[, m] = signal::filtfilt(fc$highpass, bip.montage.filt[, m]) # High pass IIR Butterwoth
}

## ----no_filtering_and_filtering_plot, include = TRUE, echo = FALSE, warning = FALSE, fig.width = 7, fig.height = 4----
ch <- 1
par(mfrow = c(2, 1), pty = "m", mai = c(0.8, 0.5, 0.5, 0.5))

# Not-filtered signal (raw signal).
plot(
  bip.montage[, ch],
  type = "l",
  panel.first = grid(),
  main = paste(colnames(bip.montage)[ch], " (raw signal, channel #1)", sep = ""),
  xlab = "sample points",
  ylab = "",
  col = "blue"
)

# Signal after filtering.
plot(
  bip.montage.filt[, ch],
  type = "l",
  panel.first = grid(),
  main = paste(colnames(bip.montage)[ch], " (filtered signal, channel #1)", sep = ""),
  xlab = "sample points",
  ylab = "",
  col = "blue"
)

par(old.par)

## ----eeg_empi_execute, include = TRUE, echo = TRUE, warning = FALSE, fig.width = 7, fig.height = 7----
# The empi.options parameter is NULL, so the EMPI program is 
# run with the parameters "-o local --gabor -i 50"

# sig <- list(bip.montage.filt, out1$sampling.rate)
# names(sig) <- c("signal", "sampling.rate")

# empi.out <- empi.execute (
#   signal = sig,
#   empi.options = NULL,
#   write.to.file = TRUE,
#   path = NULL,
#   file.name = "EEG_bipolar_filtered.db"
# )

## ----eeg_TF, include = TRUE, echo = TRUE, warning = FALSE, fig.width = 7, fig.height = 7----
# Reading a SQLite file where all the generated atom's parameters are stored.
file <- system.file("extdata", "EEG_bipolar_filtered.db", package = "MatchingPursuit")

# Generate time-frequency map based on MP atoms.
out <- empi2tf(
  db.file = file,
  channel = 1,
  mode = "sqrt",
  freq.divide = 6,
  increase.factor= 4,
  display.crosses = TRUE,
  display.atom.numbers = FALSE,
  out.mode = "plot"
)

## ----atom_params, include = TRUE, echo = TRUE, warning = FALSE----------------
file <- system.file("extdata", "sample1.db", package = "MatchingPursuit")
atom.params(file)

## ----eeg_atom_params, include = TRUE, echo = TRUE, warning = FALSE------------
file <- system.file("extdata", "EEG_bipolar_filtered.db", package = "MatchingPursuit")
ap <- atom.params(file)
head(ap, 5)
tail(ap, 5)

## ----sig_2_bin, include = TRUE, echo = TRUE, warning = FALSE------------------
file <- system.file("extdata", "sample3.csv", package = "MatchingPursuit")
out <- read.csv.signals(file)

signal.bin <- sig2bin(data = out$signal, write.to.file = FALSE)

# We have 3 channels. The first 4 time points.
head(out$signal, 4)

# The same elements of the signal in binary (floats are stored in 4 bytes).
head(signal.bin, 48)

# After decoding four sample numbers into numeric values.
# Of course we get the same values ​​as in out$signal.
readBin(signal.bin[1:4], what = "numeric", size = 4, endian = "little")
readBin(signal.bin[5:8], what = "numeric", size = 4, endian = "little")
readBin(signal.bin[41:44], what = "numeric", size = 4, endian = "little")
readBin(signal.bin[45:48], what = "numeric", size = 4, endian = "little")

## ----read_empi_db_file, include = TRUE, echo = TRUE, warning = FALSE----------
file <- system.file("extdata", "sample1.db", package = "MatchingPursuit")
out <- read.empi.db.file(file)

str(out)

## ----display_original_and_reconstrction, include = TRUE, echo = FALSE, warning = FALSE, fig.width = 7, fig.height = 4----
n.channnels <- ncol(out$original.signal)
original.signal <- out$original.signal
reconstruction <- out$reconstruction
t <- out$t
f <- out$f
len <- length(original.signal[, 1])
lab <- seq(t[1], t[len] + 1 / f, length.out = 11)

par(mfrow = c(2, 1), pty = "m", mai = c(0.8, 0.5, 0.5, 0.5))

plot(
  original.signal[,1], type = "l", col = "blue",
  main = "Original signal",
  xaxt = "n", ylab = "", xlab = "time [sec]"
)
axis(side = 1, las = 1, cex.axis = 0.9, at = seq(0, len, length.out = 11), labels = lab)

plot(
  reconstruction[,1], type = "l", col = "blue",
  main = "Rconstructed signal",
  xaxt = "n", ylab = "", xlab = "time [sec]"
)
axis(side = 1, las = 1, cex.axis = 0.9, at = seq(0, len, length.out = 11), labels = lab)

par(old.par)

## ----chirp_def, include = TRUE, echo = FALSE, warning = FALSE-----------------
n <- 1280
f <- 128
t1 = n / f
t <- seq(from = 0, to = n - 1, by = 1) / f

g1 <- gabor.fun(n, f, mean = 2, phase = 0, sigma = 1, frequency = 50, normalization = F)
g2 <- gabor.fun(n, f, mean = 8, phase = 0, sigma = 1.5, frequency = 25, normalization = F)
g3 <- gabor.fun(n, f, mean = 2, phase = 0, sigma = 2, frequency = 30, normalization = F)
g4 <- gabor.fun(n, f, mean = 7, phase = 0, sigma = 0.5, frequency = 15, normalization = F)
imp <- rep(0, n)
imp[n / 2] <- 4
sine <- 0.2 * sin(2 * pi * 3 * t)
chirp <- signal::chirp(t = t, f0 = 0, t1 = t1, f1 = 50, form = c("linear"), phase = 0)
signal <- g1$gabor + g2$gabor + g3$gabor + g4$gabor + imp  + chirp + sine

## ----chirp_plot, include = TRUE, echo = FALSE, warning = FALSE, fig.width = 7, fig.height = 7----
par(mfcol = c(8, 1), pty = "m", mai = c(0.2, 0.1, 0.2, 0.1)) 

ylim <- c(-2, 2)
plot(t, signal,   ylim = ylim, xaxt = "s", yaxt = "n", bty = "o", col = "blue",  type = "l", xlab = "", ylab = "", main = "The sum of the below functions")
plot(t, g1$gabor, ylim = ylim, xaxt = "n", yaxt = "n", bty = "o", col = "brown", type = "l", xlab = "", ylab = "",   main = "Gabor: f = 50Hz, mean = 2, sigma = 1, phase = 0")
plot(t, g2$gabor, ylim = ylim, xaxt = "n", yaxt = "n", bty = "o", col = "brown", type = "l", xlab = "", ylab = "",   main = "Gabor: f = 25Hz, mean = 8, sigma = 1.5, phase = 0")
plot(t, g3$gabor, ylim = ylim, xaxt = "n", yaxt = "n", bty = "o", col = "brown", type = "l", xlab = "", ylab = "",   main = "Gabor: f = 30Hz, mean = 2, sigma = 2, phase = 0")
plot(t, g4$gabor, ylim = ylim, xaxt = "n", yaxt = "n", bty = "o", col = "brown", type = "l", xlab = "", ylab = "",   main = "Gabor: f = 15Hz, mean = 7, sigma = 0.5, phase = 0")
plot(t, imp,      ylim = ylim, xaxt = "n", yaxt = "n", bty = "o", col = "brown", type = "l", xlab = "", ylab = "",   main = "Unit impulse")
plot(t, sine,     ylim = ylim, xaxt = "n", yaxt = "n", bty = "o", col = "brown", type = "l", xlab = "", ylab = "",   main = "Sine wave: 3 Hz")
plot(t, chirp,    ylim = ylim, xaxt = "n", yaxt = "n", bty = "o", col = "brown", type = "l", xlab = "", ylab = "",   main = "Linear-frequency chirp (0Hz - 50Hz)")

par(old.par)

## ----chirp_plot_tf, include = TRUE, echo = FALSE, warning = FALSE, fig.width = 7, fig.height = 6----
file <- system.file("extdata", "sample2.db", package = "MatchingPursuit")

out <- empi2tf(
  db.file = file,
  channel = 1,
  mode = "sqrt",
  freq.divide = 1,
  increase.factor= 4,
  display.crosses = TRUE,
  display.atom.numbers = FALSE,
  out.mode = "plot",
  plot.signals = FALSE
)


## ----Gabor_fun, include = TRUE, echo = FALSE, fig.width = 7, fig.height = 5, fig.align = 'center'----
N <- 512
fs <- 256
# normalization = T --> signal's norm = 1

par(mfrow = c(2,2), pty = "m", mai = c(0.4, 0.4, 0.3, 0.2))

n <- 4

sigmas <- c(0.5, 0.2, 0.8, 0.5)
frequencies <- c(14, 8, 4, 1)
phases <- c(0, 1, 1.5, -2)
means = c(0.5, 0.8, 1, 1.5)

s <- rep(0, N)

for (i in c(1:n)) {
  sigma <- sigmas[i]
  frequency <- frequencies[i]
  phase <- phases[i]
  mean <- means[i]
  main <- latex2exp::TeX(paste(
    "$\\mu=$", means[i], ", ", 
    "$\\sigma=$", sigmas[i], ", ",
    "$\\f=$", frequencies[i], ", ",
    "$\\phi=$", phases[i], 
    sep = ""
    )
  )
  
  gb <- gabor.fun(N, fs, mean, phase, sigma, frequency, normalization = F)

  plot(
    gb$t, 
    gb$gauss, type="l", ylim = c(-1, 1), col = "red", 
    xlab= "", ylab= "",
    xaxt = "t", yaxt = "t", bty = "or",
    cex.axis = 1, lwd = 2,
    main = main)
  lines(gb$t, gb$cosinus, type = "l", col = "grey")
  lines(gb$t, gb$gabor, col="blue", lwd = 2)
  
  s <- s + gb$gabor
}

par(old.par)

## ----restore_par, include = FALSE---------------------------------------------
par(old.par)

