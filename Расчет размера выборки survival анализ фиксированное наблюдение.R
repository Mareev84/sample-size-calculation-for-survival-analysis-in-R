# =============================================================
# Расчёт размера выборки для анализа выживаемости (лог-ранк, Шёнфельд)
# фиксированное время наблюдения за пациентами, без учёта потерь к наблюдению
# =============================================================

## ===== ПАРАМЕТРЫ ПО УМОЛЧАНИЮ =====
alpha  <- 0.05
power  <- 0.80
ratio  <- 1       # соотношение n_t / n_c (1:1); n_t — терапия, n_c — контроль
p_seq  <- seq(0.10, 0.50, by = 0.05)  # базовая доля событий в контроле за фиксированный период наблюдения
effects <- c(0.10, 0.15, 0.20, 0.25)  # эффекты: RRR (относительное снижение риска), в долях

## ===== ПРЕДПОЛОЖЕНИЯ (важно) =====
# - Фиксированное наблюдение: каждый пациент наблюдается фиксированное время (например, 6 мес) или до события.
# - Потери к наблюдению не учитываются (идеальный сценарий).
# - Экспоненциальная модель + пропорциональные риски (PH).
# - Связь долей событий за фиксированный период и HR:
#     HR = 1 - RRR
#     Если S_c = 1 - p_control, то S_t = S_c^HR => p_treated = 1 - (1 - p_control)^HR
# - Число событий по Шёнфельду переводим в n_total через ожидаемую долю событий в обеих группах.

## ===== ФУНКЦИИ =====
# Доля событий в терапии при заданных p_control и HR:
p_treated_from_pc_hr <- function(p_control, HR) 1 - (1 - p_control)^HR

# Требуемое число событий по Шёнфельду (лог-ранк, PH):
required_events_schoenfeld <- function(HR, alpha = 0.05, power = 0.80) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta  <- qnorm(power)
  ((z_alpha + z_beta)^2) / (log(HR)^2)
}

# Перевод числа событий в общий размер выборки (n_total):
# r = n_t / n_c (по умолчанию 1 для 1:1). В знаменателе — ожидаемая доля событий в популяции исследования.
required_n_total <- function(p_control, HR, alpha = 0.05, power = 0.80, r = 1) {
  D  <- required_events_schoenfeld(HR, alpha, power)
  pt <- p_treated_from_pc_hr(p_control, HR)
  n  <- D * (1 + r) / (p_control + r * pt)  # это n_total (всего в исследовании)
  ceiling(n)
}

# Удобная обёртка: подаём RRR (в долях), внутри считаем HR = 1 - RRR
required_n_from_rrr <- function(p_control, RRR, alpha = 0.05, power = 0.80, r = 1) {
  required_n_total(p_control, 1 - RRR, alpha, power, r)
}

## ===== ГРАФИК 1: n vs базовая доля событий (разные эффекты) =====
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

grid1 <- do.call(rbind, lapply(effects, function(rrr){
  data.frame(
    p_control = p_seq,
    RRR = rrr,
    n   = sapply(p_seq, function(pc) required_n_from_rrr(pc, rrr, alpha, power, ratio))
  )
}))
grid1$RRR_lab <- paste0(round(grid1$RRR*100), "%")

p1 <- ggplot(grid1, aes(x = p_control*100, y = n, color = RRR_lab)) +
  geom_line(size = 1) +
  geom_point(size = 1.8) +
  labs(
    title = "Требуемое n vs базовая доля событий (разные эффекты)",
    x = "Доля событий в контроле за фиксированный период, %",
    y = "Требуемый общий размер выборки (n)",
    color = "Эффект (RRR, %)"
  ) +
  theme_bw()
print(p1)

## ===== ГРАФИК 2: p_control фиксирован, отметить n при RRR 10/15/20/25% =====
p_control_fixed <- 0.30  # интерпретируется как доля событий за выбранный фиксированный период
df2 <- data.frame(
  RRR = effects,
  n   = sapply(effects, function(rrr) required_n_from_rrr(p_control_fixed, rrr, alpha, power, ratio))
)
df2$RRR_pct <- df2$RRR*100

p2 <- ggplot(df2, aes(x = RRR_pct, y = n)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = n), df2, linetype = "dotted", linewidth = 0.5) +
  geom_text(aes(label = n), vjust = -0.6, size = 3.8) +
  scale_x_continuous(breaks = df2$RRR_pct) +
  labs(
    title = "n при p_control=30% (α=0.05, мощность=80%, 1:1; фиксированное наблюдение — одинаковый период)",
    x = "Ожидаемое относительное снижение риска, % (RRR)",
    y = "Требуемый общий размер выборки (n)"
  ) +
  theme_bw()
print(p2)

## ===== ТАБЛИЦА ДЛЯ ПРОТОКОЛА =====
# Это n_total (всего в РКИ). Для 1:1 на группу ~ n_total/2.
tbl <- transform(df2,
                 `RRR, %` = RRR_pct,
                 `n (всего)` = n,
                 `≈ n/группу` = ceiling(n/2))
tbl <- tbl[, c("RRR, %", "n (всего)", "≈ n/группу")]
print(tbl, row.names = FALSE)


