library(ggplot2)
library(data.table)
library(hsdr)
library(patchwork)

in_prisma <- data.table::fread("WIP/data/Arborea_24082019_TOArad_PRISMA.dat")
in_6s     <- data.table::fread("WIP/data/Arborea_24082019_TOA6S_rad.dat")
wvl_6s    <- data.table::fread("WIP/data/wavelengthSE.dat")

in_6s$wvl <- wvl_6s[1:999]

n_esu <- dim(in_6s)[2]

in_6s <- in_6s[complete.cases(in_6s),]
# Divide by 100 the prisma data
in_prisma[,2:6] <- in_prisma[,2:(n_esu+1)] / 100
names(in_prisma)[2:6] <- paste("prisma", "esu", 1:n_esu, sep = "_")
names(in_6s) <- paste("6s", "esu", 1:n_esu, sep = "_")

# bind and convert to long format
rads <- t(as.matrix(in_6s[,1:n_esu]))
wls <-  in_6s$wvl
responses <- data.frame(center = in_prisma$wav, fwhm = 10)

newSpeclib <- speclib(rads, wls)
spectral_data_resampled <- spectralResampling(newSpeclib, responses, rm.NA= T, response_function = T)

resamp_data <- t(spectral_data_resampled@spectra@spectra_ma)
resamp_data <- as.data.frame(cbind(spectral_data_resampled@wavelength, resamp_data))
names(resamp_data) <- c("wvl", paste("6s", "esu", 1:n_esu, sep = "_"))

# join the datasets
specdata <- cbind(in_prisma, resamp_data[,2:(n_esu + 1)])
colprisma = paste("prisma_esu_", 1:n_esu, sep = "")
col6s = paste("6s_esu_", 1:n_esu, sep = "")
DT_melted <- data.table::melt(specdata, measure = list(colprisma, col6s),
                              variable.name = "ESU", value.name = c("PRISMA", "6S"))
DT_melted2 <- DT_melted
DT_melted2$diff     <- DT_melted$PRISMA - DT_melted$`6S`
DT_melted2$diffperc <- 100 * (DT_melted$PRISMA - DT_melted$`6S`) / DT_melted$PRISMA

p1 <- ggplot(DT_melted2) +
    geom_point(aes(x = `6S`, y = PRISMA, color = wav, pch = ESU)) +
    xlim(0,130) + ylim(0,130) +
    theme_light() +
    xlab("Simulated At Sensor Radiance") +
    ylab("PRISMA measured At Sensor Radiance") +
    scale_colour_viridis_c() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")

p2 <- ggplot(DT_melted2) +
    geom_line(aes(x = wav, y = diff, color = ESU, group = ESU)) +
    theme_light() +
    xlab("Wavelength") +
    ylab("Difference between PRISMA and Simulated At Sensor Radiance")

p3 <-ggplot(DT_melted2) +
    geom_line(aes(x = wav, y = diffperc, color = ESU, group = ESU)) +
    theme_light() +
    xlab("Wavelength") +
    ylab("% Difference between PRISMA and Simulated At Sensor Radiance") +
    ylim(-200, 100)

p4 <- ggplot(DT_melted2) +
    geom_line(aes(x = `wav`, y = PRISMA, color = ESU)) +
    theme_light() +
    xlab("Wavelength") +
    ylab("PRISMA At Sensor Radiance")

p5 <- ggplot(DT_melted2) +
    geom_line(aes(x = `wav`, y = `6S`, color = ESU)) +
    theme_light() +
    xlab("Wavelength") +
    ylab("6S At Sensor Radiance")

p4 + p5

DT_melted3 <- data.table::melt(DT_melted, id.vars = c("wav", "ESU"),
                               value.name = "Radiance")
p6 <- ggplot(DT_melted3) +
    geom_line(aes(x = `wav`, y = Radiance, color = variable)) +
    theme_light() +
    facet_wrap(~ESU, labeller = label_both) +
    xlab("Wavelength") +
    ylab("At Sensor Radiance")

save(p1,p2,p3,p4,p5,p6, file = "WIP/data/Arborea_24082019_plots.RData")
